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
