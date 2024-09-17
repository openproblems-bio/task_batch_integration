# Task Template

This repo is a template to create a new task for the OpenProblems v2. This repo contains several example files and components that can be used when updated with the task info.

> [!WARNING] 
> This README will be overwritten when performing the `create_task_readme` script.

## Create a repository from this template

> [!IMPORTANT] 
> Before creating a new repository, make sure you are part of the OpenProblems task team. This will be done when you create an issue for the task and you get the go ahead to create the task.
> For more information on how to create a new task, check out the [Create a new task](https://openproblems.bio/documentation/create_task/) documentation.

The instructions below will guide you through creating a new repository from this template ([creating-a-repository-from-a-template](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-repository-from-a-template#creating-a-repository-from-a-template)).


* Click the "Use this template" button on the top right of the repository.
* Use the Owner dropdown menu to select the `openproblems-bio` account.
* Type a name for your repository (task_...), and a description.
* Set the repository visibility to public.
* Click "Create repository from template".

## Clone the repository

To clone the repository with the submodule files, you can use the following command:

```bash
git clone --recursive git@github.com:openproblems-bio/<repo_name>.git
```
>[!NOTE]
> If somehow there are no files visible in the submodule after cloning using the above command. Check the instructions [here](common/README.md).

## What to do next

Check out the [instructions](https://github.com/openproblems-bio/common_resources/blob/main/INSTRUCTIONS.md) for more information on how to update the example files and components. These instructions also contain information on how to build out the task and basic commands.

For more information on the OpenProblems v2, check out the [documentation](https://openproblems.bio/documentation/).