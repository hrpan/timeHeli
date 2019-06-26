const int slice_types = 4;

const char *slice_vars[slice_types] = { "ep", "ed", "dist", "dt" };
const int slice_max = 7;
const int slices[slice_types] = { 7, 7, 7, 7 };

const double slice_range[slice_types][2] = {
	{ 0.7, 12 },
	{ 6, 12 },
	{ 0, 5000 },
	{ 1, 200 }
};
