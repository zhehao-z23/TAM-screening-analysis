---
name: git-server-version-control
description: Server-side Git and GitHub version-control workflow for research/code projects. Use when Codex needs to initialize or clean up a Git repository on a remote server, prepare a GitHub upload, configure HTTPS or SSH authentication, avoid committing raw data or generated results, handle an existing empty GitHub repository, or document safe branch/commit/push practices.
---

# Git Server Version Control

## Core Rule

Protect data first, then make history clean. Never upload raw sequencing data, signed download links, credentials, private SSH keys, large generated outputs, or intermediate binary files unless the user explicitly asks and the repository is designed for that.

## Project Audit

Before initializing or pushing, inspect the project:

```bash
pwd
find . -maxdepth 3 -type f | sort | head -200
git status --short --branch 2>/dev/null || true
git remote -v 2>/dev/null || true
```

Identify what belongs in Git:

- Commit: scripts, notebooks intended as templates, configs, docs, small examples, environment files, `.gitignore`.
- Ignore: FASTQ/FASTA/BAM/SAM, tar/zip archives, pickle/counts, processed data, figures/results, logs, credentials, signed URLs.
- If unsure whether a file is sensitive or huge, leave it out and ask.

## .gitignore Baseline

For sequencing/data-analysis projects, add or verify rules like:

```gitignore
# raw sequencing / archives
*.fq
*.fastq
*.fq.gz
*.fastq.gz
*.fa
*.fasta
*.bam
*.sam
*.tar
*.tar.gz
*.zip

# intermediate/generated data
*.p
*.pkl
*.pickle
processed_data/*
results/*
figures/*
data/*
!processed_data/.gitkeep
!results/.gitkeep
!figures/.gitkeep
!data/.gitkeep

# credentials and local state
.env
*.pem
*.key
id_ed25519*
known_hosts
__pycache__/
.ipynb_checkpoints/
```

Validate ignore behavior before commit:

```bash
git check-ignore -v data/test.fastq.gz results/test.png processed_data/N5_counts.p raw.tar.gz || true
```

## Initialize Repository

Use `main` as the stable branch. Use date-description branches for changes only when needed, e.g. `20260430-add-n5-qc`.

```bash
git init
git branch -m main
git config user.name "<github-user>"
git config user.email "<github-email>"
git status --short
```

Make a focused first commit:

```bash
git add .
git status --short
git commit -m "Initialize analysis workflow"
```

If files staged include raw data or results, stop and fix `.gitignore` before committing.

## Choose Authentication

Prefer SSH on servers. HTTPS often fails in non-interactive Codex/remote sessions with:

```text
fatal: could not read Username for 'https://github.com': No such device or address
```

Use HTTPS only if the user will run the push interactively with a GitHub Personal Access Token.

## SSH Key Setup

First inspect existing SSH state without printing private keys:

```bash
ls -ld ~ ~/.ssh 2>/dev/null || true
ls -l ~/.ssh 2>/dev/null || true
test -r ~/.ssh/id_ed25519 && echo private_key_readable || echo private_key_missing_or_not_readable
test -r ~/.ssh/id_ed25519.pub && echo public_key_readable || echo public_key_missing_or_not_readable
```

If `~/.ssh` is writable, generate a normal key:

```bash
ssh-keygen -t ed25519 -C "<github-email>"
cat ~/.ssh/id_ed25519.pub
```

Tell the user to add the public key to:

```text
GitHub -> Settings -> SSH and GPG keys -> New SSH key
```

If `~/.ssh` is not writable or belongs to another owner, create a dedicated key outside the repository and configure only this repo to use it:

```bash
mkdir -p /work/<user>/.github_ssh
chmod 700 /work/<user>/.github_ssh
ssh-keygen -t ed25519 -C "<github-email>" -f /work/<user>/.github_ssh/id_ed25519_github -N ""
chmod 600 /work/<user>/.github_ssh/id_ed25519_github
chmod 644 /work/<user>/.github_ssh/id_ed25519_github.pub
cat /work/<user>/.github_ssh/id_ed25519_github.pub
```

Never print or commit the private key. The private key is the file without `.pub`.

Configure the repository to use that key and a writable `known_hosts` file:

```bash
git remote set-url origin git@github.com:<owner>/<repo>.git
git config core.sshCommand "ssh -i /work/<user>/.github_ssh/id_ed25519_github -o IdentitiesOnly=yes -o UserKnownHostsFile=/work/<user>/.github_ssh/known_hosts -o StrictHostKeyChecking=accept-new"
```

Test authentication:

```bash
ssh -i /work/<user>/.github_ssh/id_ed25519_github \
  -o IdentitiesOnly=yes \
  -o UserKnownHostsFile=/work/<user>/.github_ssh/known_hosts \
  -o StrictHostKeyChecking=accept-new \
  -T git@github.com
```

Success looks like:

```text
Hi <user>! You've successfully authenticated, but GitHub does not provide shell access.
```

## Add Remote And Push

For SSH:

```bash
git remote add origin git@github.com:<owner>/<repo>.git
git push -u origin main
```

For HTTPS:

```bash
git remote add origin https://github.com/<owner>/<repo>.git
git push -u origin main
```

If remote already exists:

```bash
git remote set-url origin <url>
git remote -v
```

## Remote Already Has Content

If push is rejected with `fetch first`, do not force push by default. Inspect first:

```bash
git fetch origin main
git log --oneline --graph --decorate --all --max-count=20
```

If the remote contains real work, merge or rebase safely:

```bash
git pull --rebase origin main --allow-unrelated-histories
git push -u origin main
```

If the user confirms the remote is only an empty GitHub shell repo, use the safer force form:

```bash
git fetch origin main
git push --force-with-lease -u origin main
```

Use `--force-with-lease`, not bare `--force`, so Git refuses to overwrite if the remote changed after the last fetch.

## Validation Checklist

After push, verify:

```bash
git status --short --branch
git log --oneline --decorate --max-count=3
git remote -v
```

Expected clean state:

```text
## main...origin/main
```

Open or report the GitHub URL. Summarize the commit hash and what was included.

## Branch And Version Naming

Use:

```text
main
YYYYMMDD-short-description
```

Examples:

```text
main
20260430-init-project
20260430-add-fixed6-analysis
20260430-add-n5-qc
```

Use tags for stable analysis versions:

```bash
git tag v2026.04.30-fixed6
git push origin v2026.04.30-fixed6
```

## Research Project Notes

For analysis projects, include a `README.md` with:

- Purpose of the workflow.
- Expected input folder structure.
- Exact commands to run extraction, analysis, QC, and plotting.
- Which files are generated and ignored.
- Environment setup with `requirements.txt` or `environment.yml`.
- Data policy: raw sequencing data and generated results are not committed.

Prefer committing reusable scripts and notebook templates over one-off generated notebooks. If notebooks are committed, clear large outputs unless the output is essential documentation.
