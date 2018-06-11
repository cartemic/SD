# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2018-05-10
### Added
- _version.py
- CHANGELOG.md
- setup.cfg
- setup.py

### Changed
- README.md
- LICENSE


## [0.1.a5] - 2018-05-10
### Changed
- version number to make twine happy

## [0.1.a6] - 2018-05-17
### Changed
- added state output to detonations.cj_speed()

## [0.1.1a0-4] - 2018-05-21
### Changed
- removed garbage from __init__.py
- added package directory to setup.py
- tried more things to get conda install to work correctly

## [0.1.1a5] - 2018-05-25
### Fixed
- Equilibrium shock reflection state calculation was fixed and checked. Results are within 0.5% of the original function at ~50% of the time. Frozen shock reflection state calculation has a major bug, and has been commented out until it can be fixed.

## [0.1.1a6] - 2018-06-11
### Fixed
- Fixed problem where cj speed calculation would loop indefinitely if R^2 didn't converge above 0.9999
