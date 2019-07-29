import csv
from typing import List, Dict

def write_config_to_csv(settings: List[Dict[str, str]], filename: str) -> None:
    settings = sorted(settings, key = lambda d: d['Name'])
    with open(filename, 'x', newline='') as f:
        dw = csv.DictWriter(f, settings[0].keys())
        dw.writeheader()
        dw.writerows(settings)

def read_config_from_csv(filename: str) -> List[Dict[str, str]]:
    with open(filename, 'r') as csvfile:
        rows = [r for r in csv.reader(csvfile, delimiter=',', quotechar="\"")]
    header_row, param_rows = rows[0], rows[1:]
    return [dict(zip(header_row, r)) for r in param_rows]
