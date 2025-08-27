import os

def confirm_overwrite(file_path):
    filename = os.path.basename(file_path)
    if os.path.exists(file_path):
        response = input(f"File {filename} already exists. Do you want to overwrite it? (y/n): ").strip().lower()
        if response != 'y':
            print("Skipping...")
            return False
    return True