from codecs import open
from os import path
import sys
from setuptools import setup

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'sd2', '_version.py')) as version_file:
    exec(version_file.read())

with open(path.join(here, 'README.md')) as readme_file:
    readme = readme_file.read()

with open(path.join(here, 'CHANGELOG.md')) as changelog_file:
    changelog = changelog_file.read()

desc = readme + '\n\n' + changelog
try:
    # noinspection PyPackageRequirements
    import pypandoc
    long_description = pypandoc.convert_text(desc, 'rst', format='md')
    with open(path.join(here, 'README.rst'), 'w') as rst_readme:
        rst_readme.write(long_description)
except (ImportError, OSError, IOError):
    long_description = desc

install_requires = [
    'numpy',
    'scipy',
    'cantera'
]

tests_require = [
    'pytest',
    'pytest-cov',
]

needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
setup_requires = ['pytest-runner'] if needs_pytest else []

# noinspection PyUnresolvedReferences
setup(
    name='sd2',
    version=__version__,
    description='A SDToolbox fork for python',
    longdescription=desc,
    author='Mick Carter',
    author_email='cartemic@oregonstate.edu',
    url='https://github.com/cartemic/SD2',
    license='BSD-3-Clause',
    python_requires='>=3.6.*',
    packages=['sd2', 'sd2.tests'],
    package_dir={'sd2': 'sd2'},
    install_requires=install_requires,
    tests_require=tests_require,
    setup_requires=setup_requires
    )
