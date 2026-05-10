# pboost 0.4.3

* remove `lasso-glm`, `lasso-rq`, `lasso-sar`.

* add `penalty-glm`, `penalty-ncvreg`, `penalty-rq`, `penalty-sar`.

* add `pcoxph()`, `fcoxph()`



# pboost 0.4.2

* fix the bug of the new design version of `pboost()` (fortunately, it has not been officially released yet).

* polish `frs()`


# pboost 0.4.1

* update `lasso_rq()`.



# pboost 0.4.0

* add `lasso_sar()`, `lasso_rq()`, `lasso_glm()`.


# pboost 0.3.8

* add `coef.sarpboost()`, `set_rook_matrix()`


# pboost 0.3.7

* polish


# pboost 0.3.6

* simplify unnecessary `provided_args`



# pboost 0.3.5

* add `frq()`


# pboost 0.3.4

* add `fbetareg()`
* fix bug in`frs()`


# pboost 0.3.3

Reformulate according the new interface `pboost()`:
* `pbetareg()`
* `pglm()`
* `plm()`
* `prq()`

Newly add:
* `psar()`, `fsar()`


# pboost 0.3.2

* reformulate `pboost()`.


# pboost 0.3.1

* add `flm()`, `fglm()`.


# pboost 0.3.0

* add forward regression selection framework (`frs.R`).


# pboost 0.2.2

* update documentation.



# pboost 0.2.1

* polish the documentation.


# pboost 0.2.0

* remove commented codes in examples.

* wrap the lengthy examples in `\donttest{}`.

* eval argument `data` at first in wrapped models.



# pboost 0.1.1

* polish the help documentation.

* Refine the `DESCRIPTION` file to ensure it conforms to the required specifications



# pboost 0.1.0

* Initial CRAN submission.
