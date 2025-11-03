from pathlib import Path

def main():
    frip_files = Path("/zata/zippy/ramirezc/atac-smk/results/frip_all/0.05").glob("*.txt")
    FRIP_SCORES = []
    for file in frip_files:
        try:
            with open(file, "r") as f:
                content = f.read()
                score = float(content.split(":")[1].strip())
                FRIP_SCORES.append(score)
        except Exception as e:
            print(f"Error reading {file}: {e}")
    print(max(FRIP_SCORES))

if __name__ == "__main__":
    main()
