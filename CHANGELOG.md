# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]

### Changed

- Changed field names for `SimpleMedia` to `permittivity` and `permeability`, and added
a new two-argument constructor that derives `c`.

### Removed

- Removed `NULL_CHARGE` and `NULL_CURRENT` since new keyword-argument constructor for
`RadiationSource` allows null sources to be left unspecified.
- Converted exported `t′` function to non-exported internal `_t′` utility function.
- Removed specialty Unicode accessors for structs defined in this package.


## [0.5.x]

### Added

- Implemented GitHub CI Actions for unit testing, CompatHelper, TagBot, and a Documenter preview cleanup script.
