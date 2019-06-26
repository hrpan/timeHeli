const int slice_types = 4;

const char *slice_vars[slice_types] = { "ep", "ed", "dist", "dt" };
const int slice_max = 10;
const int slices[slice_types] = { slice_max, slice_max, slice_max, slice_max };

const double slice_range[slice_types][2] = {
	{ 1.5, 12 },
	{ 1.9, 2.7 },
	{ 0, 500 },
	{ 1, 400 }
};

