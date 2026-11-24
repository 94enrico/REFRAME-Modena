"""Build the two private review extracts from the historical panel and the updated workbook."""

from argparse import ArgumentParser
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
KEY = ["patient_id", "access_date"]
DIRECT = [
    "smoking",
    "NRTI",
    "NNRTI",
    "CD4_absolute__cells/uL",
    "CD4/CD8__ratio",
    "HIV_viral_load__copies/mL",
    "menopausa",
    "INSTI_GEN_1",
    "INSTI_GEN_2",
]
MASS = [
    "dexa_total_weight__kg",
    "dexa_total_bone_mass__g",
    "dexa_total_lean_mass__g",
    "dexa_total_fat_mass__g",
    "VAT_dexa",
]
EXAM = ["VAT_dexa", "dexa_total_fat_mass__g", "dexa_total_weight__kg"]
PANEL = [
    "sex",
    "study_end_age",
    "Nadir_CD4__cells/uL",
    "INSTI",
    "TAF",
    "IP",
    "HOMA",
    "SPPB",
    "comorbidity_COPD",
    "comorbidity_DM",
    "comorbidity_HTN",
    "comorbidity_Kidn",
    "comorbidity_CVD",
    "comorbidity_Dislipidemia",
    "WHtR",
    "waist_circumference__cm",
    "BMI",
    "phys_act_bin",
    "VAT_dexa",
    "cap_med",
    "e_med",
    "FEGMIL",
    "FIB4",
    "hiv_duration_yrs",
    "excess_adiposity",
    "insulin_resistance",
    "sppb_low",
    "bmi_cat",
]


def linked_exam(panel, raw, date, columns):
    """Use the latest assessment up to 30 days after the visit, as in the original R code."""
    examinations = raw[[*KEY, date, *columns]].copy()
    examinations["measured"] = examinations[date].fillna(examinations.access_date)
    examinations = examinations.dropna(subset=["measured"]).drop_duplicates(
        ["patient_id", "measured"]
    )
    visits = panel[KEY].assign(until=panel.access_date + pd.Timedelta(days=30)).reset_index()
    linked = pd.merge_asof(
        visits.sort_values("until"),
        examinations.drop(columns="access_date").sort_values("measured"),
        by="patient_id",
        left_on="until",
        right_on="measured",
        direction="backward",
    )
    return linked.set_index("index").sort_index()[[date, *columns]]


def make_extracts(panel, raw):
    # The original panel stores the last visit's covariates on a separate death row.
    panel = panel.loc[panel.sex.isin(["F", "M"])].reset_index(drop=True)
    death = panel.state.eq(4)
    date = panel.access_date.where(~death, panel.death_date)
    visits = panel[["patient_id", *PANEL]].copy()
    visits.insert(1, "date", date)
    visits.insert(2, "age_days", (date - panel.birth_date).dt.days)
    visits.insert(3, "adult", date.ge(panel.birth_date + pd.DateOffset(years=18)))
    visits.insert(4, "death", death)

    direct = raw.dropna(subset=["birth_date"])[[*KEY, *DIRECT]].drop_duplicates()
    extra = panel[KEY].merge(
        direct, on=KEY, how="left", sort=False, validate="many_to_one", indicator=True
    )
    if not extra["_merge"].eq("both").all():
        raise ValueError("A historical visit has no corresponding workbook record.")
    for name in DIRECT[:-2]:
        visits[name] = extra[name]
    for g in (1, 2):
        visits[f"g{g}"] = extra[f"INSTI_GEN_{g}"].map({"SI": 1, "NO": 0})
    dexa = linked_exam(panel, raw, "DEXA_date", EXAM)
    for name in dexa:
        visits[f"exam_{name}"] = dexa[name]
    hips = linked_exam(panel, raw, "antropo_date", ["waist/hips__ratio"])["waist/hips__ratio"]
    visits["hip_ratio"] = hips.where(hips.gt(0) & ~death)
    return visits, raw[["patient_id", "DEXA_date", *MASS]]


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("panel", type=Path, help="Historical REFRAME_data.csv")
    parser.add_argument(
        "workbook", type=Path, help="Updated final 4 workbook with both INSTI generations"
    )
    args = parser.parse_args()
    panel = pd.read_csv(
        args.panel,
        parse_dates=["birth_date", "access_date", "death_date"],
        float_precision="round_trip",
    )
    fields = [*KEY, "birth_date", "antropo_date", "waist/hips__ratio", "DEXA_date", *MASS, *DIRECT]
    raw = pd.read_excel(args.workbook, usecols=fields, na_values=["NULL"], engine="openpyxl")
    # One access date is stored as text in the updated workbook.
    for column in ["access_date", "birth_date", "antropo_date", "DEXA_date"]:
        raw[column] = pd.to_datetime(raw[column], format="mixed", dayfirst=True).dt.normalize()
    visits, dexa = make_extracts(panel, raw)
    destination = ROOT / "Dati"
    destination.mkdir(exist_ok=True)
    visits.to_csv(destination / "visite.csv.gz", index=False)
    dexa.to_csv(destination / "dexa.csv.gz", index=False)
