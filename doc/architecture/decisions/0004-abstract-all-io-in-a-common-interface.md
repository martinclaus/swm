# 4. Abstract all IO in a common interface

Date: 2023-01-27

## Status

Accepted

Amends [3. Apply Dependency Inversion Principle (DIP)](0003-apply-dependency-inversion-principle.md)

Amended by [6. IO shall provide a reader/writer API](0006-io-shall-provide-a-reader-writer-api.md)

## Context

Until now, IO is limited to reading and writing netCDF files to disk. In the future, alternative ways of doing IO may be required.

## Decision

To be able to easily add alternative IO methods (disk only for now), the `Io` component will be made abstract and extending types will provide concrete implementations.

## Consequences

Since the rest of the code base is depending on the abstract IO type, the dependency to a IO implementation is inverted.
This allows to plug-in alternative implementations for IO.

An important type related to (disk) IO is the `IoHandle` which represents a variable provided or received by the Io component to read and write data.
Since `IoHandle` construction may require implementation specific details (e.g. filename and variable name for disk IO, host identifier and variable name for network IO), the arguments for the constructor must be provided in a polymorphic way and their validity must be checked at runtime by the concrete IO implementation.
