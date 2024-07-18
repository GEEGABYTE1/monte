<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/RocketPy-Team/RocketPy/master/docs/static/RocketPy_Logo_white.png">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/RocketPy-Team/RocketPy/master/docs/static/RocketPy_Logo_black.png">
  <img alt="RocketPy Logo" src="https://raw.githubusercontent.com/RocketPy-Team/RocketPy/master/docs/static/RocketPy_Logo_black.png">
</picture>

<br>

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/RocketPy-Team/rocketpy/blob/master/docs/notebooks/getting_started_colab.ipynb)
[![PyPI](https://img.shields.io/pypi/v/rocketpy?color=g)](https://pypi.org/project/rocketpy/)
[![Documentation Status](https://readthedocs.org/projects/rocketpyalpha/badge/?version=latest)](https://docs.rocketpy.org/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/RocketPy-Team/RocketPy/graph/badge.svg?token=Ecc3bsHFeP)](https://codecov.io/gh/RocketPy-Team/RocketPy)
[![Contributors](https://img.shields.io/github/contributors/RocketPy-Team/rocketpy)](https://github.com/RocketPy-Team/RocketPy/graphs/contributors)
[![Chat on Discord](https://img.shields.io/discord/765037887016140840?logo=discord)](https://discord.gg/b6xYnNh)
[![Sponsor RocketPy](https://img.shields.io/static/v1?label=Sponsor&message=%E2%9D%A4&logo=GitHub&color=%23fe8e86)](https://github.com/sponsors/RocketPy-Team)
[![Instagram](https://img.shields.io/badge/Instagram-E4405F?style=flat&logo=instagram&logoColor=white)](https://www.instagram.com/rocketpyteam)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-0077B5?style=flat&logo=linkedin&logoColor=white)](https://www.linkedin.com/company/rocketpy)
[![DOI](https://img.shields.io/badge/DOI-10.1061%2F%28ASCE%29AS.1943--5525.0001331-blue.svg)](http://dx.doi.org/10.1061/%28ASCE%29AS.1943-5525.0001331)

<br>

# monte 

<div align="center">
  <img src="./content/image1.png" width="600" height="400" />
</div>

A local GUI for rocket-trajectory calculated with monte-carlo sims through RocketPy.

## running the GUI

1.Navigate to
```sh
/docs/notebooks/monte_carlo_analysis/script.py
```

2. In the terminal and within the same directory run the following command:
```sh
streamlit run script.py
```

<div align="center">
  <img src="./content/image2.png" width="600" height="400" />
</div>

The current GUI runs on the files natively given in the `notebooks/monte_carlo_analysis` directory, but in the future, will be dynamic where the user can upload data files from their 6DOF automatically. Currently, the user may need to drag the file manually into the directory.

The user is encouraged to read the `monte_carlo` tutorial from RocketPy which can be found here: https://github.com/RocketPy-Team/RocketPy/blob/master/docs/notebooks/monte_carlo_analysis/monte_carlo_analysis.ipynb 

## next steps
- dynamic file importing right from GUI
- Dynamic Environment Input instead of current static environment of Valetudo.
- Database of various "parts" of rocket so the user doesn't have to type.
- Analogous to the one directly on top, create rocket presets
- Add more computations.








