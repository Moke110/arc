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

The full native Tahoe build completed under the same 10 GiB virtual cap in
1:38:26 at 5,588,400 KiB peak RSS. WSL reported 24,612,504 KiB total memory
(23.47 GiB). The builder aggregated 95,624,334 cells across 3,388 immutable
shards into 67,018 bulk profiles on 62,710 native genes. The HDF5 contains
4,328,037,534 bytes and SHA256
`e67bfff5299878933f506f395b4bca0f4e9085af3f3cc96b6b2c456e3332cf80`.

The 128-row parity audit matched all 67,018 ordered observations, their full
deposited profile metadata and every count on the fixed 18,164-gene projection
against the sealed historical raw release. It used 415,232 KiB peak RSS and
took 1:17.46. The current native corpus verification additionally replays all
62,710 native columns directly from source.

Independent Tahoe replay and candidate publication completed in 2:06:59 under
the same cap with 7,393,328 KiB peak RSS. Source replay compared all 14 plates
across every native column, covering 67,018 profiles and 95,624,334 source cells.
The corpus publisher verified its native count payloads, coverage mappings and
source observations before replacing the release directory. The compact
[acceptance receipt](acceptance.json) records the resulting manifest and payload
digests.
