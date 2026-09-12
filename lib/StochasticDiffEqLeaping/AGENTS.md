# Jump representations

Keep mass-action and general `RegularJump` solver paths separate. Specialize
`StochasticDiffEqCore.jump_noise_data` to initialize structured reaction data;
do not normalize `MassActionJump` into `RegularJump` in the problem constructor.
Use JumpProcesses' public mass-action operations for propensities, stoichiometric
updates, and implicit drift evaluation.
