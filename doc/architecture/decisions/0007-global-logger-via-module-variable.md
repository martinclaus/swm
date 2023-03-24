# 7. Global logger via module variable

Date: 2023-03-24

## Status

Accepted

## Context

Currently, the logger object is injected into every component and passed along to all types whenever logging is required. This creates a lot of boilerplate and overhead.

## Decision

Use a module variable for logging. If a logger is initialized, it will deallocate the current value and assign itself to this variable. 

## Consequences

Instead of keeping and passing down references to a logger component, the globally accessible logger is readily available everywhere. However, only a single logging strategy can be used throughout the application, which is fine for now. Also, at least one logger needs to be created to not have undefined behaviour.  