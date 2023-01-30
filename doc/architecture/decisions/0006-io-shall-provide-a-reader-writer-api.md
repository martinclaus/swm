# 6. IO shall provide a reader/writer API

Date: 2023-01-30

## Status

Accepted

Amends [4. Abstract all IO in a common interface](0004-abstract-all-io-in-a-common-interface.md)

## Context

Having separate handel objects for reading and writing the data makes it more clear which IO is intended when creating a handle.  

## Decision

Two abstract types for reading and writing are created in `io_module`.
Those types will provide all functionality required to read or write data, respectively.

## Consequences

It is possible to inject the IO component into the reader/writer object which makes it easier to pass around those objects.
This is a change at the very abstract level which requires considerable refactoring.
