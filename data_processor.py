#data_processor.py
import pandas as pd


def read_data(file_path):
    return pd.read_csv(file_path)


def get_conversion_rules():
    return {
        'GE-primary_roots': {'Period': 'Germination & Establishment_phase(0-45days)', 'Organ': 'GE-primary_roots'},
        'GE-bud': {'Period': 'Germination & Establishment_phase(0-45days)', 'Organ': 'GE-bud'},
        'GE-tillers': {'Period': 'Germination & Establishment_phase(0-45days)', 'Organ': 'GE-tillers'},
        'GE-root_system': {'Period': 'Germination & Establishment_phase(0-45days)', 'Organ': 'GE-root_system'},
        'GE-stem': {'Period': 'Germination & Establishment_phase(0-45days)', 'Organ': 'GE-stem'},
        'GE-leaf_blades': {'Period': 'Germination & Establishment_phase(0-45days)', 'Organ': 'GE-leaf_blades'},
        'TC-root-2': {'Period': 'Tillering_phase & Canopy_development(45-120days)', 'Organ': 'TC-root'},
        'TC-stem-2': {'Period': 'Tillering_phase & Canopy_development(45-120days)', 'Organ': 'TC-stem'},
        'TC-leaf-2': {'Period': 'Tillering_phase & Canopy_development(45-120days)', 'Organ': 'TC-leaf'},
        'Ge-root-3': {'Period': 'Grand_growth & elongation_phase(120-250days)', 'Organ': 'Ge-root'},
        'Ge-stem-3': {'Period': 'Grand_growth & elongation_phase(120-250days)', 'Organ': 'Ge-stem'},
        'Ge-leaf-3': {'Period': 'Grand_growth & elongation_phase(120-250days)', 'Organ': 'Ge-leaf'},
        'MP-root-4': {'Period': 'Maturing & Ripening_phase(250-360days)', 'Organ': 'MP-root'},
        'MP-stem-4': {'Period': 'Maturing & Ripening_phase(250-360days)', 'Organ': 'MP-stem'},
        'MP-leaf-4': {'Period': 'Maturing & Ripening_phase(250-360days)', 'Organ': 'MP-leaf'},
        'MP-flower-4': {'Period': 'Maturing & Ripening_phase(250-360days)', 'Organ': 'MP-flower'},
        'MP-bud-4': {'Period': 'Maturing & Ripening_phase(250-360days)', 'Organ': 'MP-bud'},
        'MP-leaf_sheath': {'Period': 'Maturing & Ripening_phase(250-360days)', 'Organ': 'MP-leaf_sheath'},
        'MP-tillers': {'Period': 'Maturing & Ripening_phase(250-360days)', 'Organ': 'MP-tillers'},
    }

def transform_data(data_a, rules_dict):
    transformed_data = pd.DataFrame(columns=['Period', 'Organ', 'Gene', 'Expression', 'Functions'])

    print(transformed_data)

    for index, row in data_a.iterrows():
        for period_organ, rule in rules_dict.items():
            period = rule['Period']
            organ = rule['Organ']
            if period_organ in row:
                expression = row[period_organ]
            else:
                continue
            functions = row['Function']
            gene = row['Gene_Name']

            new_row = pd.DataFrame({
                'Period': [period],
                'Organ': [organ],
                'Gene': [gene],
                'Expression': [expression],
                'Functions': [functions],
                })

            transformed_data = pd.concat([transformed_data, new_row], ignore_index=True)
    return transformed_data
