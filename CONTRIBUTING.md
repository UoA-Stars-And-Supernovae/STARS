Contributing to Aotearoa STARS
==============================

Welcome to Aotearoa STARS!

Before Contributing
-------------------

1. Discuss your changes either via email or GitHub issues (preferred).
2. Build consensus that this change will be good and acceptable.
3. Scope your changes -- identify what and where will be changed.
4. Branch your changes -- we tag versions associated with individual BPASS releases, and individual projects. It is therefore of paramount importance that you appropriately branch and prep your changes.
5. Write meaningful commit messages. "Fix bugs" is bad. "Fix the CEE prescription" or "Fix #163" is better.

Getting ready to pull request
-----------------------------

0. Test your code against an observation.
1. Write a *comprehensive* PR. Describe what you've changed and why.
2. Request review by someone with write access, for both code compliance and correctness.
3. When ready, you will be asked if your PR is ready for merge, and then we will merge.

Programming: The Zen of Aotearoa STARS
--------------------------------------

IMPLICIT is the work of the devil. All programs should use IMPLICIT NONE.

We are no longer in the days where variable names were six characters.
Use full words where possible, or known abbreviations where needed.
    - DoCEE is acceptable.
    - DoCommonEnvelopeEvolution is bad.
    - DoComEnvEvo is worst.

Spaces v. Tabs. Use 4 spaces. Bind TAB to input 4 spaces if you need to.

Whitespace. One-third whitespace, one-third code, one-third comments.