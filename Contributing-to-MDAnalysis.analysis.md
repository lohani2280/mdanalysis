*DRAFT* (2015-12-08) — orbeckst

MDAnalysis is open source and welcomes contributions from the community. In fact, one of the strengths of MDAnalysis is its openness to contributors from a wide range of scientific areas, which leads to it being useful across different disciplines within a coherent framework and user interface. This is particularly reflected in the `MDAnalysis.analysis` sub-package, which contains a number of unique solutions to analysis problems such as polymer physics, water dynamics or lipid membrane physics.

In order to maintain a uniform high standard across the whole MDAnalysis library, we will introduce basic standards for code in `MDAnalysis.analysis`. We will work with contributors to achieve these standards and we will work on applying these standards to existing code.

### Standards for code in `MDAnalysis.analysis` ###

1. New module must have a **maintainer** (i.e. someone who is willing to fix bugs and review pull requests)
   * The maintainer's *name*, *email* and *github username* must be contained in the code and in the docs.
   * The date of addition must be documented (through `.. versionadded:: 1.0` reST tags and `:year: 2016` at the top of the docs).
   * The maintainer has to updated contact details (through a PR or issue).
   * The maintainer has to find a successor if they cannot maintain code anymore. The successor is then added to the docs as above.
2. Code must be fully **documented**.
   * Include all relevant *academic citations*.
   * Follow the [[Style Guide]] for doc strings.
3. Code must follow the **[[Style Guide]]**.
4. The module must have **[[Unit Tests|UnitTests]]**:
   * tests key functionality
   * comes with necessary data files
   * coverage ≥85% (as measured by coveralls)


