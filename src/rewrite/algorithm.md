```
normalize(R); normalize(C)              // cheap, deterministic, local rules only
repeat:
    P = match_prefix(R, C)              // frontier peeling from the front
    S = match_suffix(R, C)              // frontier peeling from the back
    R_mid, C_mid = what's left
    rewrite the miter  C_mid · R_mid†   // rewrite only here
    if miter is empty → proven equivalent, stop
    C_mid = rewritten C part            // cancellations may expose new matches
until P and S stop growing
```