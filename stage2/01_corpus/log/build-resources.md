# Native count resource measurements

On 2026-11-23 (Asia/Shanghai), the first Tahoe native build reached the end of
plate 1 under a 10 GiB virtual memory cap. Conversion of the full selected
plate to CSR attempted to allocate another 383 MiB of coordinate data and hit
the cap. It exited in 9m54s at 6,522,976 KiB peak RSS before sealing a checkpoint.

The native conversion now counts nonzeros and packs CSR in 64-row blocks.
A 5,000 by 62,710 count fixture with 104,520,000 nonzero values fit inside the
same cap at 3,525,684 KiB peak RSS and completed in 7.40 seconds. The projected
Tahoe builder and historical payloads retain their original bytes.

The bounded tenv and other cardiac source verification run took 46.40 seconds
at 587,072 KiB peak RSS. Its tenv replay independently checked 4,510 wells and
19,613 native genes. The raw HDF5 and probe mapping both reproduced their
recorded historical digests.

The bounded retry sealed plate 1 in 344.29 seconds as a 260 MiB checkpoint.
A 128-row streaming comparison checked all 4,643 observation identities,
5,299,930 source cells and 18,164 common columns against the accepted historical
Tahoe raw counts. Every projected count matched. The comparison took 12.14
seconds and 408,240 KiB peak RSS. The checkpoint SHA256 is
`621a8e727be14e8bb82841f59033c1cbc03ce56717bee00ee8e4bba3702a416a`.
The historical Tahoe HDF5 SHA256 is
`a830a34a11759d060fb31697df9281ad6f7667b8051616e1bfea4fbeea64ea15`.
This check covers one plate's common projection. The native release gate
independently replays all native genes on all plates.
