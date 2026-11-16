# Kv2.1 interface-protected masking pilot

This pilot is a controlled masking comparison to the registered 73-residue
Kv2.1 run. It preserves the ten extra VSD-edge positions from the 99-residue
cohort, but leaves experimentally observed cross-subunit contact positions
unmasked. The protected positions were defined from the union of heavy-atom
contacts within 5 Å in 8SD3 and 8SDA, before inspecting the pilot output.

The candidate mask is:

```text
256-260,276-280,288-294,296-297,300-301,303-304,307,309-311,
316-317,319-326,373,384,416-417
```

This is a 41-residue mask. It retains the VSD-focused additions `256-260`
and `276-280` from the 99-residue cohort while protecting the selectivity
filter, S4–S5/pore interface, and most of the S6 inter-subunit contact
network. The exact position list and contact-selection rationale are part of
the pilot provenance; this is not a replacement for the registered 73-residue
cohort until the two runs have been compared.

All three conditions use the existing Kv2.1 ColabFold settings: 100 seeds,
random seed 33, saved recycles, and 10 recycles. The input A3Ms are generated
from the unchanged vanilla A3Ms with the validated `targetedMasking` CLI and
`--channel-mask-percent 1.0` over the candidate positions.
