"""
This module contains the run function to run an executable.

When adding this file as a resource to a component, you can build the
target directory as a Python package. This will allow you to import the
run function in a script to run the executable.

Example:
  from target.executable.my_namespace.my_component import run as my_component

  out = my_component(
    input="path/to/my_input.h5ad",
    arg1="value1",
    output="output.h5ad",
    publish_dir="path/to/my_output_dir"
  )

  print(out.__dict__)
"""

import re
import os
import yaml
import glob
import tempfile
import subprocess
from pathlib import Path
from typing import List, Optional, Union
from collections import namedtuple
from dataclasses import dataclass

@dataclass
class Argument:
  type: str
  name: str
  clean_name: str
  label: Optional[str] = None
  summary: Optional[str] = None
  description: Optional[str] = None
  info: Optional[dict] = None
  default: Optional[any] = None
  example: Optional[any] = None
  required: bool = False
  direction: str = "input"
  multiple: bool = False
  multiple_sep: Optional[str] = None
  must_exist: bool = False
  create_parent: bool = False

def __read_config(dir: str) -> dict:
  """
  Read and return a Viash config file

  This function will read the config file and return the parsed config.
  It will also resolve any paths in the config's build info.

  Args:
    dir: The directory to read the config from.

  Returns:
    The parsed Viash config.
  """

  # read config
  config_path = os.path.join(dir, ".config.vsh.yaml")
  with open(config_path, "r") as f:
    config = yaml.safe_load(f)
  executable_path = os.path.join(dir, config["name"])

  # resolve paths
  build_info = config.get("build_info", {})
  build_info["output"] = dir
  build_info["config"] = config_path
  build_info["executable"] = executable_path
  config["build_info"] = build_info

  # process arguments
  config["all_arguments"] = [
    __process_argument(arg)
    for arg_group in config["argument_groups"]
    for arg in arg_group["arguments"]
  ]

  return config


def __process_argument(argument: dict) -> dict:
  """
  Process a single argument and return a cleaned version.

  This function will remove the leading dashes from the argument name and
  add a "clean_name" field. It will also add a default value if it is not
  present.

  Args:
    argument: The argument to process.

  Returns:
    The processed argument.
  """
  argument["clean_name"] = argument["name"].removeprefix("-").removeprefix("-")

  default = argument.get("default")
  type = argument["type"]
  direction = argument.get("direction", "input")

  if default is None:
    # boolean_true and boolean_false are special cases
    if type == "boolean_true":
      default = False
    elif type == "boolean_false":
      default = True
    # if it's an output file, we use the example to generate a default
    elif type == "file" and direction == "output":
      example = argument.get("example", "output.txt")

      if isinstance(example, list):
        example = example[0]

      example_basename = os.path.basename(example)
      example_ext = os.path.splitext(example_basename)[1]
      is_mult = "_*" if argument.get("multiple", False) else ""
      default = f"{argument['clean_name']}{is_mult}{example_ext}"
    argument["default"] = default

  return Argument(**argument)


# dynamically create function which has the same signature as the arguments in arguments
def __check_type_of_argument_value(value: any, argument: Argument) -> any:
  """
  Check the type of an argument and convert it if necessary.

  This function will check if the value is of the correct type and convert it
  if necessary. It will raise an error if the type is incorrect.

  Args:
    argument: The argument to check the type of.
    value: The value to check.

  Returns:
    The converted value.
  """

  if argument.type == "file" and argument.direction == "input":
    if isinstance(value, str):
      value = Path(value)
    if not isinstance(value, Path):
      raise ValueError(f"Argument {argument.clean_name} must be of type Path")
    if argument.must_exist and not value.exists():
      raise FileNotFoundError(f"File {value} not found")
  if argument.type == "file" and argument.direction == "output":
    if isinstance(value, Path):
      value = str(value)
    if not isinstance(value, str):
      raise ValueError(f"Output argument {argument.clean_name} must be a template of type str")
    if argument.multiple and not "*" in value:
      raise ValueError(f"Output argument {argument.clean_name} must contain '*' for multiple files")
  if argument.type == "integer":
    if not isinstance(value, int):
      raise ValueError(f"Argument {argument.clean_name} must be of type int")
  if argument.type == "long":
    if not isinstance(value, int):
      raise ValueError(f"Argument {argument.clean_name} must be of type int")
  if argument.type == "double":
    if not isinstance(value, float):
      raise ValueError(f"Argument {argument.clean_name} must be of type float")
  if argument.type == "boolean" or argument.type == "boolean_true" or argument.type == "boolean_false":
    if not isinstance(value, bool):
      raise ValueError(f"Argument {argument.clean_name} must be of type bool")
  return value


def __get_argument_value(argument: Argument, arg_values: dict) -> Optional[List[any]]:
  """
  Check if an argument is present and of the correct type.

  This function will check if the argument is present in the argument values
  and if it is of the correct type. It will raise an error if the argument is
  missing or if the type is incorrect.

  Args:
    argument: The argument to check.
    arg_values: The argument values to check.

  Returns:
    The converted value or None if the argument is not present.
  """
  if argument.required and argument.clean_name not in arg_values:
    raise ValueError(f"Missing required argument: {argument.clean_name}")

  value = arg_values.get(argument.clean_name)

  # if the value is None, we can skip the rest of the checks
  if value is None:
    if argument.required:
      raise ValueError(f"Argument {argument.clean_name} is required")
    return None

  # turn all values into lists for easier processing
  if not isinstance(value, list):
    value = [value]

  # check if the argument is a multiple argument
  if not argument.multiple and len(value) > 1:
    raise ValueError(f"Argument {argument.clean_name} does not accept multiple values")

  if argument.multiple:
    # todo: fix when https://github.com/viash-io/viash/issues/749 is implemented
    if argument.type == "boolean_true" or argument.type == "boolean_false":
      raise ValueError(f"Argument {argument.clean_name} does not accept multiple values")

  # check types do some type conversion
  return [__check_type_of_argument_value(v, argument) for v in value]

def __argument_value_to_flag(value: Optional[List[any]], argument: Argument) -> List[str]:
  """
  Convert an argument value to a command line flag.

  This function will convert an argument value to a command line flag. It will
  return a list of strings that can be appended to the command line.

  Args:
    value: The value to convert.
    argument: The argument to convert.

  Returns:
    The list of strings to append to the command line.
  """
  if value is None:
    # todo: update this when https://github.com/viash-io/viash/issues/705 is implemented
    return []

  if argument.type == "boolean_true":
    return [argument.name] if value[0] else []

  if argument.type == "boolean_false":
    return [] if value[0] else [argument.name]

  out = []

  for val in value:
    # todo: also update this when https://github.com/viash-io/viash/issues/705 is implemented
    if val:
      out.append(argument.name)
      out.append(str(val))

  return out

def __process_output_value(output_path: str, argument: Argument) -> Optional[Union[Path, List[Path]]]:
  """
  Process an output value and return the path(s).

  This function will process an output value and return the path(s). It will
  raise an error if the path does not exist and is required.

  Args:
    output_path: The output path to process.
    argument: The argument that the output path belongs to.

  Returns:
    The path(s) or None if the path does not exist.
  """
  if argument.multiple:
    value = glob.glob(output_path)
  else:
    value = [output_path] if os.path.exists(output_path) else []

  if argument.required and argument.must_exist and len(value) == 0:
    s = "s" if argument.multiple else ""
    raise FileNotFoundError(f"Output{s} for argument {argument.clean_name} not found at '{output_path}'")

  # convert types
  paths = [Path(v) for v in value]

  if argument.multiple:
    return paths
  else:
    return paths[0] if len(paths) > 0 else None


def __run(arg_values: dict, config: dict, publish_dir: Path, verbose: bool = False) -> dict:
  """
  Run the executable with the provided arguments.

  This function will run the executable with the provided arguments. It will
  return a dictionary with the resolved output paths.

  Args:
    arg_values: The argument values to use.
    config: The parsed Viash config.
    publish_dir: The directory to publish the output files to.

  Returns:
    The resolved output paths.
  """
  cmd = [config["build_info"]["executable"]]
  outputs = {}

  # create publish directory if it doesn't exist
  if not os.path.exists(publish_dir):
    os.makedirs(publish_dir)

  for argument in config["all_arguments"]:
    value = __get_argument_value(argument, arg_values)
    if verbose:
      print(f"Argument {argument.clean_name}: {value}", flush=True)

    if value is not None:
      # use publish dir to write output files to based on the provided template
      if argument.type == "file" and argument.direction == "output":
        value = [os.path.join(publish_dir, v) for v in value]
        outputs[argument.clean_name] = value[0] # should always be a single file

      cmd.extend(__argument_value_to_flag(value, argument))

  # run the command
  if verbose:
    print(f"Running command: {' '.join(cmd)}", flush=True)
  out = subprocess.run(cmd)

  assert out.returncode == 0, f"Command failed with return code {out.returncode}"

  resolved_outputs = {
    argument.clean_name: __process_output_value(outputs[argument.clean_name], argument)
    for argument in config["all_arguments"]
    if argument.type == "file" and argument.direction == "output"
  }

  return resolved_outputs


def __generate_function_argument_signature(argument: Argument, direction_input: bool = True) -> str:
  """
  Generate a function argument signature.

  This function will generate a function argument signature based on the
  provided argument. It will return the signature as a string.

  Args:
    argument: The argument to generate the signature for.
    direction_input: Whether we're generating the signature for an input argument or as a return value.

  Returns:
    The generated signature.
  """
  base = None
  if argument.type == "string":
    base = "str"
  elif argument.type == "file":
    if argument.direction == "input" or not direction_input:
      base = "Path"
    else:
      base = "str"
  elif argument.type == "integer":
    base = "int"
  elif argument.type == "long":
    base = "int"
  elif argument.type == "double":
    base = "float"
  elif argument.type == "boolean_true":
    base = "bool"
  elif argument.type == "boolean_false":
    base = "bool"
  elif argument.type == "boolean":
    base = "bool"
  else:
    raise ValueError(f"Unknown argument type: {argument.type}")

  if argument.multiple:
    base = f"List[{base}]"

  if not argument.required:
    base = f"Optional[{base}]"

  signature = f"{argument.clean_name}: {base}"

  if direction_input:
    if argument.default is not None or not argument.required:
      signature = f"{signature} = {repr(argument.default)}"
  else:
    if not argument.required:
      signature = f"{signature} = None"

  return signature

def __strip_margin(text: str, symbol: str = "\\|") -> str:
  """
  Strip margin from text

  This function will strip the margin from the text. It will return the text
  with the margin stripped.

  Args:
    text: The text to strip the margin from.
    symbol: The margin symbol to use.

  Returns:
    The text with the margin stripped.

  Example:
    __strip_margin(
      '''\
      |Hello
      |World
      |'''
    )
    # returns "Hello\nWorld\n"
  """

  return re.sub("(^|\n)[ \t]*" + symbol, "\\1", text)

def __generate_output_dataclass(config: dict) -> str:
  """
  Generate the output dataclass for the executable.

  This function will generate the output dataclass for the executable based on
  the provided config. It will return the dataclass as a string

  Args:
    config: The parsed Viash config.

  Returns:
    The generated output dataclass.
  """
  outputs = [
    __generate_function_argument_signature(arg, direction_input=False)
    for arg in config["all_arguments"]
    if arg.direction == "output"
  ]

  body = "\n  ".join(outputs)

  return __strip_margin(
    f"""\
    |@dataclass
    |class Output:
    |  {body}
    |"""
  )

def __generate_help(config: dict) -> str:
  cmd = [config["build_info"]["executable"], "--help"]

  out = subprocess.run(cmd, capture_output=True, text=True)

  return out.stdout

def __generate_function_code(config: dict) -> str:
  """
  Generate the function code for the executable.

  This function will generate the function code for the executable based on
  the provided config. It will return the code as a string.

  Args:
    config: The parsed Viash config.

  Returns:
    The generated function code.
  """
  # order arguments by whether they have a default and then by index
  arguments = config["all_arguments"]
  arguments = sorted(arguments, key=lambda x: (x.default is not None, arguments.index(x)))

  input_arguments = [
    arg
    for arg in arguments
    if arg.direction == "input" or arg.type == "file"
  ]

  # note that output files should also be inputs
  # since the user can pass a template for the output file
  # that needs to be generated in the publish directory

  # get the signature of inputs without defaults
  inputs_no_defaults = [
    __generate_function_argument_signature(arg, direction_input=True)
    for arg in input_arguments
    if arg.default is None and arg.required
  ]
  # get the signature of inputs with defaults
  inputs_with_defaults = [
    __generate_function_argument_signature(arg, direction_input=True)
    for arg in input_arguments
    if not (arg.default is None and arg.required)
  ]

  # signature of the combined inputs
  inputs = inputs_no_defaults + [f"publish_dir: Path"] + inputs_with_defaults

  # generate help
  help = __generate_help(config)
  help_lines = help.split("\n")
  help_block = "\n  ".join(help_lines)

  return __strip_margin(
    f"""\
    |def run({', '.join(inputs)}) -> Output:
    |  '''
    |  {help_block}
    |  '''
    |  inputs_dict = locals()
    |  publish_dir = inputs_dict.pop('publish_dir')
    |  output_dict = __run(inputs_dict, config=__config, publish_dir=publish_dir)
    |  return Output(**output_dict)
    |"""
  )

__dir = os.path.dirname(__file__)
__config = __read_config(__dir)
exec(__generate_output_dataclass(__config))
exec(__generate_function_code(__config))

# # for debugging
# __dir = "viash_components/executable/differential_expression/deseq2"
# __config = __read_config(__dir)
# print(__generate_output_dataclass(__config))
# print(__generate_function_code(__config))
