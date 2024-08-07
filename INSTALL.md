# Install notes

> This project uses Miniconda to create a python environment
> <https://docs.anaconda.com/miniconda/>

- Create conda environment

  ```
  conda env create -f environment.yml -n coseismic_vertical
  ```

- Update conda environment

  ```
  conda env update -f environment.yml
  ```

- If using PyCharm,
  [configure a pycharm project with an existing conda environment](https://docs.anaconda.com/working-with-conda/ide-tutorials/pycharm/#configuring-a-pycharm-project-with-an-existing-conda-environment)
