import setuptools

# with open("README.md", "r", encoding="utf-8") as fh:
#    long_description = fh.read()

setuptools.setup(
    name="QualityFactor",
    version="0.0.1",
    author="Ofir Milul",
    author_email="ofirmilul@gmail.com",
    description="calculating quality facor for microwave resonator usin S21 response",
    # long_description=long_description,
    # long_description_content_type="text/markdown",
    #url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)