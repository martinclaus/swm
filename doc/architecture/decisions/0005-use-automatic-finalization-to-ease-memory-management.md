# 5. Use automatic finalization to ease memory management

Date: 2023-01-27

## Status

Accepted

## Context

Derived owning allocatable data (also as pointer) must release the memory when the type gets out of scope.
Since Fortran 2003, derived types can have a final binding that can be used to clean up memory upon destruction. 

## Decision

Types containing allocatable components or pointer components which gets allocated and represent data owned by the type get a final binding to clean up memory

## Consequences

Keeping track of allocated memory and freeing it up becomes easier at high level.
However, when derived types are defined, special care must be taken to properly handle resources owned by the type.
In particular, pointer components must be analysed if the data pointed to is "owned" by the type, that is if it should be deallocated when the instance is destroyed.
