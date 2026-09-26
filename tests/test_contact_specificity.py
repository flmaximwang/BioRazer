"""Cross-model contact specificity (``report_contact_specificity``)."""
import pandas as pd
import pytest

from biorazer.structure.analysis.static.report import report_contact_specificity

# Macro labels as BioRazer-PyMOL writes them -- here they are just the row keys,
# the report treats labels as opaque strings.
PAIR_SHARED = ("///A/ALA`1/N", "///B/GLY`2/O")
PAIR_ONLY_MODEL1 = ("///A/ALA`1/CB", "///B/GLY`2/CA")


def _write(path, columns, rows):
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False)
    return path


def test_specificity_index_counts_the_models_a_contact_appears_in(tmp_path):
    # One file per model, each in the column spelling its writer produced --
    # PyMOL writes ``atom1``/``distance (...), the summary uses capitals.
    model1 = _write(
        tmp_path / "model1.csv",
        ["atom1", "atom2", "distance"],
        [[*PAIR_SHARED, 3.1], [*PAIR_ONLY_MODEL1, 3.4]],
    )
    model2 = _write(
        tmp_path / "model2.csv",
        ["Atom1", "Atom2", "Distance"],
        [[*PAIR_SHARED, 3.2]],
    )

    table = report_contact_specificity([model1, model2])

    assert list(table.columns) == [
        "atom1",
        "atom2",
        "Specificity Index",
        "model1",
        "model2",
    ]
    index_of = dict(zip(zip(table.atom1, table.atom2), table["Specificity Index"]))
    assert index_of[PAIR_SHARED] == 2
    assert index_of[PAIR_ONLY_MODEL1] == 1

    shared = table[
        (table.atom1 == PAIR_SHARED[0]) & (table.atom2 == PAIR_SHARED[1])
    ].iloc[0]
    assert (shared["model1"], shared["model2"]) == (3.1, 3.2)

    # a model that lacks the contact leaves the cell empty, not 0
    only_model1 = table[table.atom1 == PAIR_ONLY_MODEL1[0]].iloc[0]
    assert pd.isna(only_model1["model2"])


def test_table_is_written_to_the_requested_csv(tmp_path):
    model1 = _write(
        tmp_path / "only.csv", ["atom1", "atom2", "distance"], [[*PAIR_SHARED, 3.1]]
    )
    output = tmp_path / "Summary.csv"

    table = report_contact_specificity([model1], output_file=output)

    assert output.exists()
    assert pd.read_csv(output).to_dict("list") == table.to_dict("list")


def test_empty_contact_list_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="empty"):
        report_contact_specificity([])

    header_only = _write(tmp_path / "empty.csv", ["atom1", "atom2", "distance"], [])
    empty = report_contact_specificity([header_only])
    assert empty.empty
    assert list(empty.columns) == [
        "atom1",
        "atom2",
        "Specificity Index",
        "empty",
    ]


def test_model_without_contacts_keeps_its_column(tmp_path):
    model1 = _write(
        tmp_path / "has.csv", ["atom1", "atom2", "distance"], [[*PAIR_SHARED, 3.1]]
    )
    model2 = _write(tmp_path / "none.csv", ["atom1", "atom2", "distance"], [])

    table = report_contact_specificity([model1, model2])

    assert list(table.columns) == [
        "atom1",
        "atom2",
        "Specificity Index",
        "has",
        "none",
    ]
    assert table["none"].isna().all()
    assert table["Specificity Index"].tolist() == [1]
