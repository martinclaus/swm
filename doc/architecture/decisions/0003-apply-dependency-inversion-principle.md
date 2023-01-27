# 3. Apply Dependency Inversion Principle (DIP)

Date: 2023-01-16

## Status

Accepted

Amends [2. Separate system into independently deployable components](0002-separate-system-into-independently-deployable-components.md)

Amended by [4. Abstract all IO in a common interface](0004-abstract-all-io-in-a-common-interface.md)

## Context

Clean architecture implies that volatile components should depend on less volatile components, or low-level components on high-level components

## Decision

For some components the dependencies must point against the flow of control.
This includes the components:
-`state`
-`io`
-`log`

## Consequences

DIP requires a mechanism, which may be Dependency Injection. Dependencies may be injected directly as polymorphic objects within the `main` component or abstract factories for those objects may be injected.
