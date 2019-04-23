const int slice_types = 4;

const char *slice_vars[slice_types] = { "ep", "ed", "dist", "dt" };
const int slice_max = 5;
const int slices[slice_types] = { 5, 5, 5, 5 };

const double slice_range[slice_types][2] = {
	{ 0.7, 12 },
	{ 6, 12 },
	{ 0, 5000 },
	{ 1, 200 }
};

