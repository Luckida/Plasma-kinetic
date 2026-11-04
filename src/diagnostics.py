"""Stable reaction-history columns with explicit gaps when a reaction leaves the top five."""
import numpy as np


def reaction_history(ranking):
    """Постоянный столбец реакции; пустая ячейка вне top-5, включая повторный вход."""
    groups = {}
    for row in ranking:
        group = groups.setdefault((row["target"], row["direction"]),
            dict(target=row["target"], target_label=row["target_label"], direction=row["direction"], columns={}, points={}))
        point = group["points"].setdefault(row["time_s"], dict(total=row["total_cm3_s"], coverage=row["top5_coverage"], reactions={}))
        if row["rank"] > 0:
            group["columns"].setdefault(row["id"], row["reaction"])
            point["reactions"][row["id"]] = row
    histories = []
    for group in groups.values():
        # Два разных механизма с одинаковой формулой не сливаются в один столбец.
        labels = list(group["columns"].values())
        columns = {ident: label if labels.count(label) == 1 else f"{label} [{ident}]"
                   for ident, label in group["columns"].items()}
        rates, fractions = [], []
        for time, point in sorted(group["points"].items()):
            common = dict(time_s=time, top5_coverage=point["coverage"])
            rate_row = dict(common, total_cm3_s=point["total"])
            fraction_row = dict(common)
            for ident, column in columns.items():
                entry = point["reactions"].get(ident)
                rate_row[column] = entry["events_cm3_s"] if entry else ""
                fraction_row[column] = entry["fraction_of_total"] if entry else ""
            rates.append(rate_row)
            fractions.append(fraction_row)
        histories.append(dict(target=group["target"], target_label=group["target_label"],
            direction=group["direction"], columns=columns, rates=rates, fractions=fractions))
    return histories
