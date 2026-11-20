# Human motif reference for Tangermeme

The AlphaGenome motif-importance stage expects a nonredundant or clustered
human transcription-factor motif database in MEME format at:

```text
data/reference/motifs/human_clustered.meme
```

The database must be selected and downloaded before running the DVC stage. Its
source, release/version, license, download URL, checksum, clustering strategy,
and motif-name convention must be recorded here. A redundant collection should
not be substituted silently because near-duplicate PWMs fragment motif counts
and complicate comparisons among EVN, SVN, and EVS.

The small MEME files distributed with Tangermeme's test suite are fixtures and
must not be used for biological interpretation.
