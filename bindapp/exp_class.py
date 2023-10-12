class Experiment:
    def __init__(self, ref_binds: dict, exp_binds: dict, scope: int):
        self.ref_binds = ref_binds
        self.exp_binds = exp_binds
        self.scope = scope
        self.num_overlaps = 0
        self.overlaps = []
        self.overlap_full = []
        self.overlap_front = []
        self.overlap_end = []
        self.overlap_ext = []
        self.font = {'fontname':'Sans Serif'}

    def within(self, exp_bind: tuple, ref_bind: tuple):
        """Returns the overlap between exp_bind and ref_bind."""
        midpoint = (int)((int(ref_bind[0]) + int(ref_bind[1])) / 2)
        # Create the scope.
        overlap = range(
            max(int(exp_bind[0]), midpoint - self.scope),
            min(int(exp_bind[1]), midpoint + self.scope) + 1,
        )
        if len(overlap) == 0:
            return overlap
        else:
            # Adjust to -1000 to 1000 base pair range.
            overlap = [x - midpoint for x in overlap]
            return overlap

    def compare_chrom_bind(self):
        """Compares the binding sites for each chromosome."""
        # Make all the Binding Comparisons
        for chromosome in self.ref_binds.keys():
            if chromosome not in self.exp_binds:
                continue
            for ref_bind in self.ref_binds[chromosome]:
                for exp_bind in self.exp_binds[chromosome]:
                    overlap = self.within(exp_bind, ref_bind)
                    if len(overlap) == 0:
                        # No overlap in the binding sites... continue.
                        continue
                    else:
                        if exp_bind[0] >= ref_bind[0] and exp_bind[1] <= ref_bind[1]:
                            # The RNA peak is fully contained by the DNA peak
                            self.overlap_full.extend(overlap)
                            overlap_type = "Fully Contained"
                            ot = "OF"
                        elif exp_bind[1] > ref_bind[1] and exp_bind[0] <= ref_bind[1]:
                            # The RNA Peak overlaps the end of the DNA peak.
                            self.overlap_end.extend(overlap)
                            overlap_type = "Overlaps with End of Ref Peak"
                            ot = "OE"
                        elif exp_bind[0] < ref_bind[0] and exp_bind[1] >= ref_bind[0]:
                            # The RNA Peak overlaps the front of the DNA peak.
                            self.overlap_front.extend(overlap)
                            overlap_type = "Overlaps with Beg. of Ref Peak"
                            ot = "OB"
                        else:
                            self.overlap_ext.extend(overlap)
                            overlap_type = "Overlaps Outside of the Ref. Peak"
                            ot = "OX"
                        self.num_overlaps += 1
                        self.overlaps.append((ref_bind, exp_bind, overlap_type, ot))

    def get_num_overlaps(self):
        """Returns the number of overlaps."""
        return self.num_overlaps

    def get_all_overlaps(self):
        """Returns all the overlaps."""
        return (
            self.overlaps,
            self.overlap_full,
            self.overlap_front,
            self.overlap_end,
            self.overlap_ext,
        )
    
    def generate_pie


class Intersection:
    def __init__(self, ref_bed: Bed, exp_bed: Bed, scope: int):
        self.ref_bed = ref_bed
        self.exp_bed = exp_bed
        self.scope = scope
        self.overlaps = []
