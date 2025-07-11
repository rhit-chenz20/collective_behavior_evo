import os

def split_all_ind_files_by_generation(root_folder):
    for dirpath, _, filenames in os.walk(root_folder):
        for fname in filenames:
            if fname.endswith(".tsv"):
                file_path = os.path.join(dirpath, fname)

                base = os.path.splitext(fname)[0]
                output_dir = dirpath  # save in the same folder

                gen_buffers = {}

                with open(file_path, "r") as infile:
                    for line in infile:
                        fields = line.strip().split("\t")
                        if len(fields) < 2:
                            continue
                        gen = fields[0]
                        if gen not in gen_buffers:
                            gen_buffers[gen] = []
                        gen_buffers[gen].append(line)
                        
                print(f"Processing {file_path} with {len(gen_buffers)} generations.")
                for gen, lines in gen_buffers.items():
                    if type(gen) is not int:
                        next
                    out_filename = f"{base}_gen{gen}.tsv"
                    out_path = os.path.join(output_dir, out_filename)
                    with open(out_path, "w") as out_file:
                        out_file.writelines(lines)
                    print(f"Created {out_path} with {len(lines)} lines for generation {gen}.")

# Run the batch split for ../data1/ind
split_all_ind_files_by_generation("/ind")
