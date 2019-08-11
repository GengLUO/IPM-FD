Inner Produck Masking with faults detection (IPM-FD) countermeasures for AES
======

We provide a C implementation of the IPM-FD countermeasure as described in [IPM-FD](#references).
The countermeasure is implemented for the AES-128 block-cipher.

What is implemented
-------------------

* AES with the  IPM-FD countermeasure without key expansion protection
* AES with the  IPM-FD countermeasure with protected Key Expansion

Notes
----

* We have not protected the key-schedule in the 1st implementation. Therefore we assume that the block-cipher initially receives the shares of the subkeys, instead of the shares of the key. Moreover we have not implemented the refresh of the key between executions; therefore the implementation would be secure only in a restricted model in which always the same intermediate variables are probed. To get security in the full model one would need to refresh the subkeys between executions.

References
----------

[IPM-FD] Wei Cheng, Claude Carlet, Kouassi Goli, Jean-Luc Danger, and Sylvain Guilley. Detecting Faults, in Inner-Product Masking Scheme - IPM-FD: IPM with Fault Detection. In 8th International Workshop on Security Proofs for Embedded Systems Atlanta, USA, August 24, 2019.
