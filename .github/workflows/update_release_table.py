import requests
import argparse
import os
import re

# GitHub repository details
# adjust REPO to match your setup
REPO = "TravisWheelerLab/ULTRA"
API_URL = f"https://api.github.com/repos/{REPO}/releases"

def get_releases_data():
    headers = {
        'Authorization': f'token {os.getenv("GITHUB_TOKEN")}',
        'Accept': 'application/vnd.github.v3+json',
    }
    response = requests.get(API_URL, headers=headers)
    response.raise_for_status()
    return response.json()

def find_release_by_version(releases, version):
    for release in releases:
        if release["tag_name"] == version:
            return release
    return None

def extract_os_arch_from_filename(filename):
    pattern = r'[.-]([a-zA-Z]+-[a-zA-Z]+)\.([a-zA-Z0-9_]+)\.tar\.gz'
    
    match = re.search(pattern, filename)
    
    if match:
        os = match.group(1)
        arch = match.group(2)

        if os == "macos-latest":
            os = "MacOS"
        elif os == "ubuntu-latest":
            os = "Ubuntu"
        elif os == "windows-latest":
            os = "Windows"

        if arch == "x64":
            arch = "Intel/AMD 64-bit"
        elif arch == "386":
            arch = "Intel/AMD 32-bit"
        elif arch == "arm64":
            if os == "MacOS":
                arch = "M1/M2/M3 (ARM 64-bit)"
            else:
                arch = "ARM 64-bit"
        elif arch == "arm":
            arch = "ARM 32-bit"

        return os, arch
    else:
        return None, None

def generate_markdown_table(release):
    table = "### Release Assets\n"
    table += "| OS | Architecture | Link |\n"
    table += "|---------|----------|-------------|\n"
    
    for asset in release["assets"]:
        if not asset["name"].endswith(".md5"):
            # parse
            os, arch = extract_os_arch_from_filename(asset["name"])
            download_url = asset["browser_download_url"]
            table += f"| {os}  | {arch}  | [Download]({download_url}) |\n"

    # Add note about Mac binary signing restriction
    table += (
        "\nTo address the Mac binary signing restriction, use the following command: "
        "`sudo xattr -dr com.apple.quarantine <path to file>/my-binary-amd64`\n"
    )
    
    return table

def update_release_body(release_id, new_body):
    headers = {
        'Authorization': f'token {os.getenv("GITHUB_TOKEN")}',
        'Accept': 'application/vnd.github.v3+json',
    }
    update_url = f"https://api.github.com/repos/{REPO}/releases/{release_id}"
    data = {
        "body": new_body
    }
    
    response = requests.patch(update_url, headers=headers, json=data)
    response.raise_for_status()

"""
def update_readme(table):
    with open("README.md", "r") as file:
        lines = file.readlines()

    with open("README.md", "w") as file:
        inside_table = False
        for line in lines:
            if line.startswith("| Version"):
                inside_table = True
                file.write(table)
                continue
            if inside_table and line.startswith("|"):
                continue
            file.write(line)
"""

def main(target_version):
    releases = get_releases_data()
    release_to_update = find_release_by_version(releases, target_version)
    
    if release_to_update:
        release_id = release_to_update["id"]
        current_body = release_to_update.get("body", "")

        markdown_table = generate_markdown_table(release_to_update)
        print(markdown_table)

        if "### Release Assets\n" not in current_body:
            updated_body = current_body + "\n\n" + markdown_table
            update_release_body(release_id, updated_body)
            print(f"Release {target_version} updated successfully.")
    else:
        print(f"Release with version {target_version} not found.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update a specific release in GitHub.")
    parser.add_argument("version", help="The tag name of the release to update")
    args = parser.parse_args()
    
    main(args.version)
