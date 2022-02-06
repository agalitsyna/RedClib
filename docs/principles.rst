Design principles of RNA-DNA nextflow pipeline
==============================================

RNA-DNA interactions data is diverse and complex. Design of RNA-DNA workflows resulted in a set of principles that
helped to generalize over a set of RNA-DNA protocols.

- One workflow is full data processing of RNA-DNA interactions from raw reads to contact pairs, binned interactions tables and annotations.
Each RNA-DNA **experimental protocol** has its own pipeline (`redc.nf` for Red-C, `gridseq.nf` for GRID-Seq, etc.).
Typically, you don't need to modify the pipeline for you data, unles it's completely new RNA-DNA protocol.

- Each **pipeline** is built out of **modules**, ready-made blocks that can be combined in different ways to
produce the results for various modifications of RNA-DNA capture. Typically, you don't need to modify modules, because
they have adjustable parameters that cover a broad range of use cases.

- We specify parameters in separate files. One run of a workflow is typically controlled by **two input files**:
configuration file (e.g., `params-redc.yml`) and sample spreadsheet (`samplesheet-redc.csv`).
For each protocol, we have separate config and sample spreadsheet, which you can modify according to your needs.

- **Config file** has important information about the RNA/DNA fragments parsing, filtering and mapping.
It typically has standardized format, which is unsafe to change. For each required parameters we specify that in the
docstring of the config. *Please, read the config docstrings carefully!*

- **Samplesheet** has information about the actual pipeline input data: path to the files with raw FASTQ reads
(produced with certain type of RNA-DNA capture protocol) or SRA/DRA/ERA ID that allows you to automatically download the
data from public databases. There, you also specify the name of the samples, library, type of sequencing (paired-end or single-end)
and read length.

- Several types of **outputs** are typically located in the `results/` folder. These include:
table(s) with contacts statistics, table(s) with individual contacts passing the filters, cool files with binned contacts
and annotated contacts. All of these outputs are optional and can be adjusted. For example, you may specify which
statistics you want to report and control the columns of the final tables with contacts.

Pipeline design
---------------
This section is more technical and might serve you well if you want to debug your RNA-DNA pipeline
or plan to adapt it to completely new type of protocol. The pipeline is powered by nextflow DSL2. Thus, it inherits the
key principles of any nextflow DSL2 pipeline: data is propagated through the channels that connect different
processes. Each process takes one or multiple files and generates one or multiple new files.

Each process corresponds to a single module. Some of the modules are joined into subworkflows, if they are
tightly connected. For example, `redc.nf` has the oligos mapping subworkflow that takes the read table and reads and
performs search of bridge/adaptors with fast Rabin-Karp algorithm. Within this subworkflow, we need to create index files for oligos
and FASTQs, then we need to parse the output of `rklib`. These steps are rather technical and might be
compartmentalized as a subworkflow, which performs a single function: search of oligos.

- **Modules and subworkflows import.** In the beginning of workflow, the modules and subworkflows are imported.
These imports have also the definitions of the parameters of the processes.
Some of the parameters are unsafe to change and are part of the particular workflow (e.g. `redc.nf` does not
allow you to analyse extra mapping tags other than NH and XM from hisat2 output).
Other parameters can be changed and are inherited from the config (and thus these parameters can be controlled by the config).

- **Typical channel structure.** Each process operates on a channel or multi-channel input.
For the consistency, one channel is typically a `[meta, file(s)]`  list. Meta is a hashMap with information about
the files and sometimes about the settings for the process.

- **Metadata.** Meta has very important key `id` that is used by any module. Other important keys include

    - `id` represents the process in the nextflow report and serves as prefix to the generated file.
By default, it's the name of the sample and replicate (inherited from the samplesheet), but might also include the
chunk number and other properties of data. When the channel expands (for example, when creating a separate instance for each DNA/RNA fragment) a special tag
is added to the `id` of its meta (e.g., name of the fragment).

    - `original_id` is an original id of the sample that is never modified (so that initial id can be always reconstructed,
for example, when we group the chunks).

    - `single_end` for type of data in the analysis.

    - `side`: 1 for forward read and 2 for reverse read, refers to the strand of origin for the fragment, or the side to which the
filter should be applied, or where the oligo should be called.

    - `fragment`: name of the fragment (after extraction step).

There are extra parameters that are added and removed simultaneously (listed in the beginning of the pipeline).

- **Channels restructuring logic.** You may notice that the most of the code of the pipeline nextflow file
is the definition of channels and manupilations on them.
This channels restructuring logic is very important part of the pipeline, and typically consists of simple steps:

    - (1) collect the data (e..g from config file)

    - (2) connect the channels

    - (3) filter the channel

    - (4) add/remove distinctive keys to the meta or tags to the id

    - (5) convert a channel to multi-channel (wwhen a process requires separate synchronized channels as an input)

    - (6) synchronize channels

The latter step of synchronization can be achieved by one of two nextflow operations combination:
(1) mix and groupTuple; (2) combine and filter. (1) is prefered for multi-files inputs (e.g. multiple tables for merge, which number
might very depending on the protocol). (2) is best for synching the finite and well-defined number of channels (e.g.,
read table combined with filters to produce fragments).

- **Resume.** One of the great advantages of nextflow is that pipelines can be stopped at any point and then
continued from the last completed step. However, this requires all the operations within the code to be deterministic.
It means that no processes can take inputs from non-synchronized channels, the files should be sorted the same way
if they come in a list through the same channel, and many more. We pay great attention to make our
pipelines deterministic and possible to resume, however, it is not ideal yet, and sometimes running
with `-resume -cache true` results in repetition of some of the computations.

Supposedly, you might try 'lenient' mode, as we noticed that for HISAT2_ALIGN the hash is frequently determined inconsistently,
probably because the list of index files is long.

Coding principles and hints
--------------------------
Most coding principles are inherited from nextflow's DSL2 and nf-core modules.
There are few initiatives that we take to simplify the code and make it more structured and readable.

- **Inline documented code.**
Usually, it's hard to read complex Groovy operations on channels.
We take the time to document each channel transformation.

- **Separate repetitive code elements into functions.**
With nextflow, there are usually repetitive operations on channels:
    - Add extra parameters to the meta hashMap with ``update_meta``.
    - Remove extra parameters of the meta with ``removeKeys``.
    - Add tag to the process id with ``add_tag``.
    - Remove one or several tags from the processid: ``remove_tag``

- Do not use ``.get(key, default_value)`` method. If it does not find the key in
the map, it creates one with ``defaultValue``. Use ``.getOrDefault(key, default_value)`` instead.

