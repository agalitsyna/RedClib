Design principles of RNA-DNA nextflow pipeline
=====

- One run of a pipeline if a single workflow. **Input:**  file with raw FASTQ reads
(produced with certain type of RNA-DNA capture protocol).
**Output:** table with contacts statistics and (optional) contacts passing filters.

- Typical channel structure

- Meta

- Required parameters of meta and extra parameters

- tags

- mix and groupTuple vs combine and filter

- Deterministic code.
Problems encountered: 'lenient' mode for HISAT2_ALIGN because the index is typically too long.

Coding principles
------
Most coding principles are inherited from nextflow's DSL2 and nf-core modules.
There are few initiatives that I take to simplify the code and make it more structured and readable.

- **Inline documented code.**
Usually, it's hard to read complex Groovy operations on channels.
We take the time to document each channel transformation.

- **Separate repetitive code elements into functions.**
With nextflow, there are usually repetitive operations on channels:
    - Add extra parameters to the meta hashMap with ``update_meta``.
    - Remove extra parameters of the meta with ``removeKeys``.
    - Add tag to the process id with ``add_tag``.
    - Remove one or several tags from the processid: ``remove_tag``

- Do not use ``.get(key, default_value)`` method! If it does not find the key in
the map, it creates one with ``defaultValue``. Use ``.getOrDefault(key, default_value)`` instead.

