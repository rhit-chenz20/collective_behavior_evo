#!/usr/bin/env python3
import argparse
from pathlib import Path
import shutil
import sys

def unique_path(p: Path) -> Path:
    """Return a non-colliding path by appending _1, _2, ... before the suffix."""
    if not p.exists():
        return p
    stem, suffix = p.stem, p.suffix
    parent = p.parent
    i = 1
    while True:
        cand = parent / f"{stem}_{i}{suffix}"
        if not cand.exists():
            return cand
        i += 1

def main():
    parser = argparse.ArgumentParser(
        description="Recursively rename files by replacing a substring in the filename and move to a new root."
    )
    parser.add_argument("src_root", type=Path,
                        help="Source root directory (e.g., n_group_extreme/data)")
    parser.add_argument("dst_root", type=Path,
                        help="Destination root directory (e.g., n_group_translational/data)")
    parser.add_argument("--find", default="ave", help='Substring to find in filename (default: "ave")')
    parser.add_argument("--replace", default="ext", help='Replacement substring (default: "ext")')
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite existing files at destination (default: off)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would happen without moving anything")
    args = parser.parse_args()

    src_root: Path = args.src_root.resolve()
    dst_root: Path = args.dst_root.resolve()

    if not src_root.is_dir():
        print(f"ERROR: Source root does not exist or is not a directory: {src_root}", file=sys.stderr)
        sys.exit(1)

    moved = 0
    renamed = 0
    unchanged = 0

    for src_path in src_root.rglob("*"):
        if not src_path.is_file():
            continue

        rel = src_path.relative_to(src_root)
        src_name = src_path.name

        # Rename basename only
        new_name = src_name.replace(args.find, args.replace)
        if new_name != src_name:
            renamed += 1
        else:
            unchanged += 1

        dst_path = (dst_root / rel.parent / new_name)

        # Ensure destination directory exists
        dst_path.parent.mkdir(parents=True, exist_ok=True)

        # Handle collisions
        final_dst = dst_path if (args.overwrite or not dst_path.exists()) else unique_path(dst_path)

        if args.dry_run:
            action = "MOVE"
            note = ""
            if final_dst != dst_path and not args.overwrite:
                note = " (collision → renamed)"
            print(f"{action}: {src_path}  ->  {final_dst}{note}")
        else:
            # Move (also acts as rename)
            shutil.move(str(src_path), str(final_dst))
            moved += 1

    if args.dry_run:
        print("\n[DRY RUN] No files moved.")
    print(f"Files renamed (basename changed): {renamed}")
    print(f"Files unchanged (no '{args.find}' in name): {unchanged}")
    if not args.dry_run:
        print(f"Files moved: {moved}")

if __name__ == "__main__":
    main()
