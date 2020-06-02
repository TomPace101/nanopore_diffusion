import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="simproc",
    version="0.1",
    author="Tom Pace",
    author_email="TomPace101@gmail.com",
    description="Simulation processing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=["simproc"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ]
)