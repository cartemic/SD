from codecs import open
from os import path
import sys

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'sd2', '_version.py')) as version_file:
    exec(version_file.read())

with open(path.join(here, 'README.md')) as readme_file:
    readme = readme_file.read()

with open(path.join(here, 'CHANGELOG.md')) as changelog_file:
    changelog = changelog_file.read()

desc = readme + '\n\n' + changelog
try:
    import pypandoc
    long_description = pypandoc.convert_text(desc, 'rst', format='md')
    with open(path.join(here, 'README.rst'), 'w') as rst_readme:
        rst_readme.write(long_description)
except (ImportError, OSError, IOError):
    long_description = desc

install_requires = [
    'numpy',
]

tests_require = [
    'pytest',
    'pytest-cov',
]

needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
setup_requires = ['pytest-runner'] if needs_pytest else []

setup(
    name='sd2',
    version=__version__,
    description='A SDToolbox fork for python',
    author='Mick Carter',
    author_email='cartemic@oregonstate.edu',
    url='https://github.com/cartemic/SD2',
    classifiers=[
        'License :: OSI Approved :: BSD License',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
    license='BSD-3-Clause',
    python_requires='>=2.7.*, >=3.5.*',
    zip_safe=False,
    package=find_packages(),
    include_package_data=True,
)

if __name__ == "__main__":
    setup()
