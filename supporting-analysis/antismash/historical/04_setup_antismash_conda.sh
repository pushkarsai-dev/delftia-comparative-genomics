#!/bin/bash
# ============================================================================
#  antiSMASH 8.0.4 — setup via conda (no Docker)
#
#  Same antiSMASH version as the Docker image, so results are equivalent.
#
#  Note on architecture: bioconda publishes antiSMASH for linux-64, noarch
#  and osx-64 -- but NOT osx-arm64. On your M2 Pro the environment is
#  therefore pinned to osx-64 and runs under Rosetta. That is the same
#  emulation cost Docker would have had.
#
#  Usage:   bash "04_setup_antismash_conda.sh"
# ============================================================================

set -u
ENVNAME="antismash8"
ASVER="8.0.4"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

say()  { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
note() { printf '    %s\n' "$*"; }
fail() { printf '\n\033[31mERROR: %s\033[0m\n\n' "$*" >&2; exit 1; }

# ------------------------------------------------------------------ conda
say "Locating conda"
# Test that conda actually RUNS -- the presence of conda.sh is not enough.
# The Anaconda at /opt/anaconda3 has a conda launcher whose python is
# missing, so it exists but cannot start. Miniforge is checked first.
CONDA_SH=""
for c in "/opt/miniforge3" "$HOME/miniforge3" "/opt/homebrew/Caskroom/miniforge/base" \
         "$HOME/anaconda3" "$HOME/opt/anaconda3" "/opt/anaconda3" \
         "$HOME/miniconda3" "$HOME/opt/miniconda3" "/opt/miniconda3"; do
  [ -f "$c/etc/profile.d/conda.sh" ] || continue
  if [ -x "$c/bin/conda" ] && "$c/bin/conda" --version >/dev/null 2>&1; then
    CONDA_SH="$c/etc/profile.d/conda.sh"; break
  else
    note "skipping broken install: $c"
  fi
done
[ -n "$CONDA_SH" ] || fail \
"No working conda found.

 Your Anaconda at /opt/anaconda3 is broken -- its conda launcher points at
 /opt/anaconda3/bin/python, which does not exist. That is why the previous
 run failed; it was never a dependency conflict.

 Install a clean, minimal conda alongside it (no admin password, does not
 touch Anaconda):

     bash 07_install_miniforge.sh

 then re-run this script."
note "using $CONDA_SH"
# shellcheck disable=SC1090
source "$CONDA_SH"

# ---------------------------------------------------------------- rosetta
say "Checking Rosetta (needed to run x86-64 packages on Apple Silicon)"
if [ "$(uname -m)" = "arm64" ]; then
  if arch -x86_64 /usr/bin/true 2>/dev/null; then
    note "Rosetta is present"
  else
    note "Rosetta is NOT installed. Installing it now -- macOS will ask you"
    note "to agree to a licence and may ask for your password."
    softwareupdate --install-rosetta --agree-to-license || fail \
"Rosetta install failed. Run this yourself and then re-run the script:
     softwareupdate --install-rosetta"
  fi
else
  note "Intel Mac, Rosetta not needed"
fi

# ------------------------------------------------------------ environment
say "Creating the '$ENVNAME' environment with antiSMASH $ASVER"
note "pinned to osx-64; this takes several minutes"

if conda env list | awk '{print $1}' | grep -qx "$ENVNAME"; then
  note "environment '$ENVNAME' already exists -- skipping creation"
  note "(to start over: conda env remove -n $ENVNAME)"
else
  CONDA_SUBDIR=osx-64 conda create -y -n "$ENVNAME" \
      -c bioconda -c conda-forge \
      "antismash=$ASVER" || fail \
"Environment creation failed. Copy the last 20 lines of the error and
 retain them for troubleshooting the dependency conflict."
fi

conda activate "$ENVNAME" || fail "Could not activate '$ENVNAME'."
# Make the architecture pin sticky for anything installed later.
conda config --env --set subdir osx-64 2>/dev/null

command -v antismash >/dev/null 2>&1 || fail \
"antiSMASH did not install into the environment. Check the installation log."

say "antiSMASH reports:"
antismash --version

# ------------------------------------------------------------- databases
say "Downloading antiSMASH databases (~9 GB, one time)"
note "this is the slow part; leave it running"
if antismash --check-prereqs >/dev/null 2>&1; then
  note "databases already present and complete -- skipping"
else
  download-antismash-databases || fail \
"Database download failed. Re-run this script -- it resumes."
fi

say "Preparing data files and caches"
antismash --prepare-data || note "(--prepare-data reported an issue; continuing)"

# ------------------------------------------------------- record provenance
say "Recording version provenance"
{
  echo "route:            conda (bioconda)"
  echo "environment:      $ENVNAME"
  echo "antismash:        $(antismash --version 2>&1 | head -1)"
  echo "conda subdir:     osx-64 (Rosetta on $(uname -m))"
  echo "host:             macOS $(sw_vers -productVersion), $(uname -m)"
  echo "created:          $(date '+%Y-%m-%d %H:%M:%S %Z')"
} | tee "$HERE/antismash_version.txt"

say "Setup complete"
cat <<EOF

  Next step -- run the panel:
      bash "$HERE/05_run_panel40_conda.sh"

EOF
