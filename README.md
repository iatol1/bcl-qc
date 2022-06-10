### Goal

Build the minimal infrastructure needed for demultiplexing and quality assurance of data generated by Illumina sequencers in a clinical NGS lab.

### Strategy

For an NGS lab like ours with few engineers - the less infrastructure we build, the less we need to maintain, a lower chance of failure, and easier to debug. Our strategy is to outsource the most complex pieces of necessary infrastructure to mature commercial solutions with a large number of users with similar use-cases. More users ensure exhaustive bug reporting and pushes vendors to prioritize fixes and new features. At the same time, we don't want to invest in "end-to-end" solutions from vendors (e.g. [ICA](https://developer.illumina.com/news-updates/illumina-connected-analytics-productionize-your-informatics-workflows-at-scale)), and must instead carefully modularize our infrastructure, so that we can replace each module with better alternatives in the future. Modular infrastructure also gives us some autonomy to solve issues ourselves without dependence on vendor tech support.

### Scope

- Setup a linux server with redundant storage that runs scheduled tasks and hosts a QC and variant sign-off GUI.
- Develop scripts to orchestrate use of [Microsoft's Azure](https://azuremarketplace.microsoft.com/en-us/marketplace/apps/illuminainc1586452220102.public-dragen-batch) and [Illumina's Dragen](https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html) for storage/compute.
- Deploy a minimal BCL to VCF pipeline to generate run and sample-level metrics for quality assurance.

### Server

[Click here](dx-var.md) to see how we set up our primary server. These instructions are written for sysadmins or DevOps engineers, but if you are a quick study and armed with Google, you'll do fine too. Detailed and versioned instructions like these helps to quickly replace the primary server in case of failure. ::TODO:: Later, we should use configuration management tools like [Ansible](https://en.wikipedia.org/wiki/Ansible_(software)) to automate deployment of new servers for fault tolerance or load balancing.

### Quick start

```bash
git clone git@github.com:ucladx/bcl-qc.git
cd bcl-qc
pip install .
```
::TODO:: Put this in PyPI instead
