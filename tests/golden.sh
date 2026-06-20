#!/bin/bash
# Golden-file parity: run the Rust dockq-rs CLI on the upstream examples and diff against
# the reference's testdata/*.dockq (header lines containing '*' are stripped, exactly as
# upstream's run_test.sh does). Long-format cases strip the header; --short cases (no
# header) compare verbatim.
set -uo pipefail

DOCKQ="/Users/Andre.Teixeira/projects/DockQ"
BIN="/Users/Andre.Teixeira/projects/dockq-rs/target/release/dockq-rs"

if [ ! -x "$BIN" ]; then
    echo "FATAL: build the release CLI first: cargo build --release -p dockq-cli"
    exit 2
fi
cd "$DOCKQ" || exit 2

pass=0
fail=0

run_case() {
    local name="$1" golden="$2" mode="$3"
    shift 3
    local out
    out="$("$BIN" "$@" 2>/dev/null)"
    local tmp
    tmp="$(mktemp)"
    echo "$out" >"$tmp"
    local a b
    if [ "$mode" = "strip" ]; then
        a="$(grep -v '\*' "$tmp")"
        b="$(grep -v '\*' "$DOCKQ/testdata/$golden")"
    else
        a="$(cat "$tmp")"
        b="$(cat "$DOCKQ/testdata/$golden")"
    fi
    rm -f "$tmp"
    if [ "$a" = "$b" ]; then
        echo "PASS $name"
        pass=$((pass + 1))
    else
        echo "FAIL $name"
        diff <(echo "$a") <(echo "$b") | head -40
        fail=$((fail + 1))
    fi
}

run_case 1A2K            1A2K.dockq                    strip examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb
run_case 1A2K_noalign    1A2K.dockq                    strip examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --no_align
run_case dimer_dimer     dimer_dimer.dockq             exact examples/dimer_dimer.model.pdb examples/dimer_dimer.pdb --short
run_case model_native    model.dockq                   strip examples/model.pdb examples/native.pdb --allowed_mismatches 1
run_case 1EXB            1EXB.dockq                    exact examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short
run_case 1EXB_AB.BA      1EXB_AB.BA.dockq              exact examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping 'AB*:BA*'
run_case 1EXB_.ABC       1EXB_.ABC.dockq               exact examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping ':ABC'
run_case 1EXB_full       1EXB_ABCDEFGH.BADCFEHG.dockq  exact examples/1EXB_r_l_b.model.pdb examples/1EXB_r_l_b.pdb --short --mapping 'ABCDEFGH:BADCFEHG'

echo ""
echo "golden: $pass passed, $fail failed"
exit $([ "$fail" -eq 0 ] && echo 0 || echo 1)
