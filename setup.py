from setuptools import setup

setup(
    name="steel_beam_analysis",
    packages=["steel_beam_analysis"],
    description="A Python program used to analyze steel beams.",
    version="0.0.0",
    url="https://github.com/jessebluestein/steel_beam_analysis",
    author="Jesse Bluestein",
    author_email="jessebluestein@gmail.com",
    keywords=[
        "steel_beam_analysis",
        "finite_element_analysis",
        "structural_engineering",
    ],
    include_package_data=True,
    install_requires=[
        "matplotlib",
        "numpy",
        "pint",
        "jinja2",
        "pdflatex",
        "sortedcontainers",
        "pytexit",
    ],
)
