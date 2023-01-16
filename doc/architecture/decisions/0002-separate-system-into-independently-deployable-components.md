# 2. Separate system into independently deployable components

Date: 2023-01-16

## Status

Accepted

## Context

Most numerical ocean models are a mess leading to excessively high costs for developing new features.
Clean architecture builds on separating components by clear architectural boundaries where dependencies always point towards less volatility.

## Decision

System is decomposed into following components:
- `main`: wires everything together.
- `state`: in-memory application state.
- `core`: Commonly used data types and interfaces, such as `Grid`, `Variable`.
- `components`: Contains abstract types declaring common interfaces for all components, e.g. `BaseComponent` with methods `initialize` and `finalize` or `DynComponent` extending `BaseComponent` with methods `step` and `advance`.
- `app`: Application controller in which the main loop runs. It orchestrates all components.
- `swm`: Solver of the shallow water equations.
- `tracer`: Solver of the tracer equation including sources and sinks.
- `log`: Implements logging to stdout.
- `io`: Implements IO to and from the local file system.
- `diag`: Online diagnostics.

## Consequences

By creating architectural boundaries between those components, code volatility will be isolated which enhances maintainability, and increased flexibility will increase extensibility. This comes at the cost of a major refactoring. Potentially, backwards compatibility will be broken requiring a new major release.
