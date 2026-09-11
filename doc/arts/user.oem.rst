Retrieval transformations
=========================

ARTS provides built-in retrieval transformations described in
:doc:`concept.oem`.  Custom
transformations can be assigned directly to any Jacobian target from Python by
providing its three operators:

.. code-block:: python

  target = ws.jac_targets.atm[-1]
  target.transform_state = lambda t, field: A @ (t - b)
  target.inverse_state = lambda x, field: A_inv @ x + b
  target.inverse_jacobian = lambda J, x, field: J @ A_inv

Here ``field`` is the complete owning field or data object, allowing mappings
that need information beyond the target itself.  Each callable must return a
vector or matrix with the same shape as the target block.  The example supports
a general invertible affine transformation; ``A`` need not be diagonal or
orthogonal.  The same interface can express bounded and other reversible
functional transformations.

See :doc:`concept.oem` for the forward/inverse transformation definitions and
the Jacobian chain rule used by these operators.
