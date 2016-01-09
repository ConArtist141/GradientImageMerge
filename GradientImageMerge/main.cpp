#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>
#include <string>
#include <algorithm>
#include <cmath>

#define _USE_MATH_DEFINES
#include <math.h>

#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <Eigen/LU>
#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>

#include "imagewin32.h"
#include "lodepng.h"
#include "GridGraph_2D_4C.h"

#define DEFAULT_IMAGE_SOURCE_1 "mountains2.png"
#define DEFAULT_IMAGE_SOURCE_2 "mountains3.png"
#define DEFAULT_IMAGE_OUTPUT "result.png"
#define DEFAULT_STITCH_MARGIN 200
#define DEFAULT_MODE PoissonStitchSparse

// It seems that setting the gradient epsilon > 0.5f causes a strange feedback effect
// The same is true for laplace and > 0.25f
#define GRADIENT_DESCENT_ITERATIONS 20
#define GRADIENT_DESCENT_EPSILON 0.5f
#define LAPLACIAN_DESCENT_ITERATIONS 40
#define LAPLACIAN_DESCENT_EPSILON 0.25f
#define RECOVERY_GRADIENT_ITERATIONS 40
#define RECOVERY_LAPLACE_ITERATIONS 80

using namespace std;

// Program mode
enum ProgramMode
{
	SimpleStitch,
	RecoverFromGradient,
	RecoverFromLaplace,
	GradientStitch,
	LaplaceStitch,
	PoissonStitch,
	PoissonStitchSparse,
	ComputeGradient,
};

struct pairhash {
public:
	template <typename T, typename U>
	std::size_t operator()(const std::pair<T, U> &x) const
	{
		return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
	}
};

// 2-component vector
template <typename t>
struct vec2
{
	t x, y;
};

template <typename t>
inline t magnitudeSquared(vec2<t> vec)
{
	return vec.x * vec.x + vec.y * vec.y;
}

template <typename t>
inline vec2<t> difference(vec2<t> a, vec2<t> b)
{
	return vec2<t> { a.x - b.x, a.y - b.y };
}

template <typename t>
inline t distanceSquared(vec2<t> a, vec2<t> b)
{
	return magnitudeSquared(difference(a, b));
}

template <typename t1, typename t2>
inline vec2<t1> scalarMult(vec2<t1> v, t2 a)
{
	return vec2<t1> { v.x * a, v.y * a };
}

template <typename t>
inline vec2<t> vecAdd(vec2<t> a, vec2<t> b)
{
	return vec2<t> { a.x + b.x, a.y + b.y };
}

typedef vector<vector<float> > ScalarField;
typedef vector<vector<vec2<float> > > VectorField;
typedef GridGraph_2D_4C<float, float, float> Grid;

template <typename t>
void copyField(const vector<vector<t> >& src, vector<vector<t> >* dest)
{
	int width = src[0].size();
	for (int y = 0, height = src.size(); y < height; ++y)
		memcpy((*dest)[y].data(), src[y].data(), sizeof(t) * width);
}

template <typename t>
void resizeField(const int width, const int height, vector<vector<t> >* field)
{
	field->resize(height);
	for (auto& vec : *field)
		vec.resize(width);
}

void computeStitchMask(ScalarField field1, ScalarField field2, ScalarField* maskResult)
{
	assert(field1.size() == field2.size());
	assert(field1[0].size() == field2[0].size());

	auto gridWidth = field1[0].size();
	auto gridHeight = field1.size();
	auto bigNumber = 1000.0f * gridWidth * gridHeight;

	resizeField(gridWidth, gridHeight, maskResult);

	Grid grid(gridWidth, gridHeight);

	// Set terminals
	for (size_t y = 0; y < gridHeight; ++y)
	{
		grid.set_terminal_cap(grid.node_id(0, y), bigNumber, 0.0f);
		grid.set_terminal_cap(grid.node_id(gridWidth - 1, y), 0.0f, bigNumber);
	}

	// Set grid
	for (size_t y = 0; y < gridHeight; ++y)
	{
		for (size_t x = 0; x < gridWidth - 1; ++x)
		{
			auto weight = fabs(field1[y][x] - field2[y][x + 1]) + fabs(field2[y][x] - field1[y][x + 1]);
			grid.set_neighbor_cap(grid.node_id(x, y), 1, 0, weight);
			grid.set_neighbor_cap(grid.node_id(x + 1, y), -1, 0, weight);
		}
	}

	for (size_t y = 0; y < gridHeight - 1; ++y)
	{
		for (size_t x = 0; x < gridWidth; ++x)
		{
			auto weight = fabs(field1[y][x] - field2[y + 1][x]) + fabs(field2[y][x] - field1[y + 1][x]);
			grid.set_neighbor_cap(grid.node_id(x, y), 0, 1, weight);
			grid.set_neighbor_cap(grid.node_id(x, y + 1), 0, -1, weight);
		}
	}

	// Compute maxflow
	grid.compute_maxflow();

	for (size_t y = 0; y < gridHeight; ++y)
	{
		for (size_t x = 0; x < gridWidth; ++x)
		{
			auto segment = grid.get_segment(grid.node_id(x, y));
			if (segment == 0)
				(*maskResult)[y][x] = 1.0f;
			else
				(*maskResult)[y][x] = 0.0f;
		}
	}
}

void computeGradientStitchMask(VectorField field1, VectorField field2, ScalarField* maskResult)
{
	assert(field1.size() == field2.size());
	assert(field1[0].size() == field2[0].size());

	auto gridWidth = field1[0].size();
	auto gridHeight = field1.size();

	maskResult->resize(gridHeight);
	for (auto& row : *maskResult)
		row.resize(gridWidth);

	float bigNumber = 1000.0f * gridWidth * gridHeight;
	Grid grid(gridWidth, gridHeight);

	// Set terminals
	for (size_t y = 0; y < gridHeight; ++y)
	{
		grid.set_terminal_cap(grid.node_id(0, y), bigNumber, 0.0f);
		grid.set_terminal_cap(grid.node_id(gridWidth - 1, y), 0.0f, bigNumber);
	}

	// Set grid
	for (size_t y = 0; y < gridHeight; ++y)
	{
		for (size_t x = 0; x < gridWidth - 1; ++x)
		{
			auto weight = distanceSquared(field1[y][x], field2[y][x + 1]) + distanceSquared(field2[y][x], field1[y][x + 1]);
			grid.set_neighbor_cap(grid.node_id(x, y), 1, 0, weight);
			grid.set_neighbor_cap(grid.node_id(x + 1, y), -1, 0, weight);
		}
	}

	for (size_t y = 0; y < gridHeight - 1; ++y)
	{
		for (size_t x = 0; x < gridWidth; ++x)
		{
			auto weight = distanceSquared(field1[y][x], field2[y + 1][x]) + distanceSquared(field2[y][x], field1[y + 1][x]);
			grid.set_neighbor_cap(grid.node_id(x, y), 0, 1, weight);
			grid.set_neighbor_cap(grid.node_id(x, y + 1), 0, -1, weight);
		}
	}

	// Compute maxflow
	grid.compute_maxflow();

	for (size_t y = 0; y < gridHeight; ++y)
	{
		for (size_t x = 0; x < gridWidth; ++x)
		{
			auto segment = grid.get_segment(grid.node_id(x, y));
			if (segment == 0)
				(*maskResult)[y][x] = 1.0f;
			else
				(*maskResult)[y][x] = 0.0f;
		}
	}
}

template <typename t>
void mergeRows(const vector<vector<vector<t> >* >& fields, vector<vector<t> >* output)
{
	auto fieldWidth = 0;
	auto fieldHeight = (*fields[0]).size();
	auto fieldCount = fields.size();

	for (auto field : fields)
		fieldWidth += (*field)[0].size();

	output->resize(fieldHeight);
	for (auto& row : *output)
		row.resize(fieldWidth);

	for (size_t y = 0; y < fieldHeight; ++y)
	{
		int xPos = 0;
		for (auto field : fields)
		{
			auto rowSrc = (*field)[y];
			memcpy(&((*output)[y].data())[xPos], rowSrc.data(), sizeof(t) * rowSrc.size());
			xPos += rowSrc.size();
		}
	}
}

template <typename t>
void subRange(const vector<vector<t> >& input, const int xStart, const int yStart, const int width,
	const int height, vector<vector<t> >* output)
{
	output->resize(height);
	for (auto& row : *output)
		row.resize(width);

	for (int y = 0, ysource = yStart; y < height; ++y, ++ysource)
		for (int x = 0; x < width; ++x)
			memcpy((*output)[y].data(), &input[ysource].data()[xStart], sizeof(t) * width);
}

void computeLaplacian(const ScalarField& scalarField, ScalarField* output)
{
	output->resize(scalarField.size());
	for (auto& row : *output)
		row.resize(scalarField[0].size());

	// Compute laplacian
	for (int y = 0, ymax = static_cast<int>(scalarField.size()); y < ymax; ++y)
		for (int x = 0, xmax = static_cast<int>(scalarField[0].size()); x < xmax; ++x)
		{
			int sampleMinusX = std::max(0, x - 1);
			int samplePlusX = std::min(xmax - 1, x + 1);
			int sampleMinusY = std::max(0, y - 1);
			int samplePlusY = std::min(ymax - 1, y + 1);
			float ddx = (scalarField[y][samplePlusX] - 2.0f * scalarField[y][x] + scalarField[y][sampleMinusX]);
			float ddy = (scalarField[samplePlusY][x] - 2.0f * scalarField[y][x] + scalarField[sampleMinusY][x]);
			(*output)[y][x] = ddx + ddy;
		}
}

// Compute a gradient from a scalar field
void computeGradient(const ScalarField& scalarField, VectorField* output)
{
	// Correct dimensions
	output->resize(scalarField.size());
	for (auto& row : *output)
		row.resize(scalarField[0].size());

	// Compute gradient
	for (int y = 0, ymax = static_cast<int>(scalarField.size()); y < ymax; ++y)
		for (int x = 0, xmax = static_cast<int>(scalarField[0].size()); x < xmax; ++x)
		{
			int sampleMinusX = std::max(0, x - 1);
			int samplePlusX = std::min(xmax - 1, x + 1);
			int sampleMinusY = std::max(0, y - 1);
			int samplePlusY = std::min(ymax - 1, y + 1);
			(*output)[y][x].x = (scalarField[y][samplePlusX] - scalarField[y][sampleMinusX]) / 2.0f;
			(*output)[y][x].y = (scalarField[samplePlusY][x] - scalarField[sampleMinusY][x]) / 2.0f;
		}
}

// Apply a mask to two images
template <typename t>
void applyMask(const vector<vector<t> >& field1, const vector<vector<t> >& field2, const ScalarField& mask, vector<vector<t> >* output)
{
	assert(field1.size() == field2.size());
	assert(field1[0].size() == field2[0].size());

	size_t width = field1[0].size();
	size_t height = field1.size();

	output->resize(field1.size());
	for (auto& vec : *output)
		vec.resize(field1[0].size());

	for (size_t y = 0; y < height; ++y)
	{
		for (size_t x = 0; x < width; ++x)
		{
			float maskValue = mask[y][x];
			(*output)[y][x] = maskValue * field1[y][x] + (1.0f - maskValue) * field2[y][x];
		}
	}
}

// Apply a mask to two images
template < >
void applyMask<vec2<float> >(const VectorField& field1, const VectorField& field2, const ScalarField& mask, VectorField* output)
{
	assert(field1.size() == field2.size());
	assert(field1[0].size() == field2[0].size());

	size_t width = field1[0].size();
	size_t height = field1.size();

	output->resize(field1.size());
	for (auto& vec : *output)
		vec.resize(field1[0].size());

	for (size_t y = 0; y < height; ++y)
	{
		for (size_t x = 0; x < width; ++x)
		{
			float maskValue = mask[y][x];
			(*output)[y][x] = vecAdd(scalarMult(field1[y][x], maskValue), scalarMult(field2[y][x], 1.0f - maskValue));
		}
	}
}

// Convert gradient data to image data
void convertGradientToImageData(const VectorField& grad, vector<unsigned char>* output)
{
	output->resize(4 * grad.size() * grad[0].size());
	int index = 0;
	for (auto& row : grad)
	{
		for (auto val : row)
		{
			(*output)[index++] = ((unsigned char)(std::fmax(0.0, val.x) * 255.0f));
			(*output)[index++] = ((unsigned char)(std::fmax(0.0, val.y) * 255.0f));
			(*output)[index++] = 0;
			(*output)[index++] = 255;
		}
	}
}

// Convert an 8-bit image array to a matrix of floats
void convertImageDataToFloatMatrix(const vector<unsigned char>& image,
	const int width, const int height, ScalarField* output)
{
	for (int y = 0; y < height; ++y)
	{
		output->push_back(vector<float>(width));
		vector<float>& row = (*output)[output->size() - 1];

		for (int x = 0; x < width; ++x)
			row[x] = (float)image[4 * (x + y * width)] / 255.0f;
	}
}

// Convert a matrix of float to an 8-bit image array
void convertFloatMatrixToImageData(const ScalarField& matrix,
	vector<unsigned char>* output)
{
	output->resize(4 * matrix[0].size() * matrix.size());
	int index = 0;
	for (const vector<float>& vec : matrix)
		for (float val : vec)
		{
			(*output)[index++] = ((unsigned char)(val * 255.0f));
			(*output)[index++] = ((unsigned char)(val * 255.0f));
			(*output)[index++] = ((unsigned char)(val * 255.0f));
			(*output)[index++] = 255;
		}
}

unsigned int floatMatrixFromPNG(const string& filename, ScalarField* output)
{
	// Open the raw PNG data
	vector<unsigned char> imageData;
	unsigned int imageWidth;
	unsigned int imageHeight;
	unsigned int error = lodepng::decode(imageData, imageWidth, imageHeight, filename);

	if (error)
		return error;

	// Convert data to a float matrix
	convertImageDataToFloatMatrix(imageData, imageWidth, imageHeight, output);
	return 0;
}

unsigned int saveFloatMatrixToPNG(const string& filename, const ScalarField& data)
{
	vector<unsigned char> outputImageData;
	convertFloatMatrixToImageData(data, &outputImageData);

	// Save the resulting image
	return lodepng::encode(filename, outputImageData, static_cast<unsigned int>(data[0].size()),
		static_cast<unsigned int>(data.size()));
}

// Solve using standard gradient descent
void gradientDescent(const ScalarField& initialGuess, const VectorField& vectorField,
	const int iterations, const float epsilon, const bool bVerticalBoundaryConstraints,
	const bool bHorizontalBoundaryConstraints, ScalarField* output)
{
	assert(initialGuess.size() == vectorField.size());
	assert(initialGuess[0].size() == vectorField[0].size());

	int width = initialGuess[0].size();
	int height = initialGuess.size();

	ScalarField tempField1;
	ScalarField tempField2;

	// Copy over initial geuss
	resizeField(width, height, output);
	resizeField(width, height, &tempField1);
	resizeField(width, height, &tempField2);

	copyField(initialGuess, &tempField1);
	copyField(initialGuess, &tempField2);

	ScalarField* readField = &tempField1;
	ScalarField* writeField = &tempField2;

	int yStart = (bHorizontalBoundaryConstraints ? 1 : 0);
	int yEnd = (bHorizontalBoundaryConstraints ? height - 1 : height);
	int xStart = (bVerticalBoundaryConstraints ? 1 : 0);
	int xEnd = (bVerticalBoundaryConstraints ? width - 1 : width);

	for (int i = 0; i < iterations; ++i)
	{
		cout << "Descent iteration " << i + 1 << "..." << endl;

		for (int y = yStart; y < yEnd; ++y)
		{
			// The left and right boundaries are left untouched
			for (int x = xStart; x < xEnd; ++x)
			{
				auto xm = max(x - 1, 0);
				auto xmm = max(x - 2, 0);
				auto xp = min(x + 1, width - 1);
				auto xpp = min(x + 2, width - 1);
				auto ym = max(y - 1, 0);
				auto ymm = max(y - 2, 0);
				auto yp = min(y + 1, height - 1);
				auto ypp = min(y + 2, height - 1);
				auto df = 0.5f * (4.0f * (*readField)[y][x] - (*readField)[y][xpp] - (*readField)[y][xmm] - (*readField)[ypp][x] - (*readField)[ymm][x])
					+ vectorField[y][xp].x - vectorField[y][xm].x + vectorField[yp][x].y - vectorField[ym][x].y;
				(*writeField)[y][x] = fmax(fmin((*readField)[y][x] - epsilon * df, 1.0f), 0.0f);
			}
		}

		// Swap the fields
		auto temp = writeField;
		writeField = readField;
		readField = temp;
	}

	copyField(*readField, output);
}

// Solve using laplacian descent
void laplacianDescent(const ScalarField& initialGuess, const ScalarField& laplacian,
	const int iterations, const float epsilon, const bool bVerticalBoundaryConstraints,
	const bool bHorizontalBoundaryConstraints, ScalarField* output)
{
	assert(initialGuess.size() == laplacian.size());
	assert(initialGuess[0].size() == laplacian[0].size());

	int width = initialGuess[0].size();
	int height = initialGuess.size();

	ScalarField tempField1;
	ScalarField tempField2;

	// Copy over initial geuss
	resizeField(width, height, output);
	resizeField(width, height, &tempField1);
	resizeField(width, height, &tempField2);

	copyField(initialGuess, &tempField1);
	copyField(initialGuess, &tempField2);

	ScalarField* readField = &tempField1;
	ScalarField* writeField = &tempField2;

	int yStart = (bHorizontalBoundaryConstraints ? 1 : 0);
	int yEnd = (bHorizontalBoundaryConstraints ? height - 1 : height);
	int xStart = (bVerticalBoundaryConstraints ? 1 : 0);
	int xEnd = (bVerticalBoundaryConstraints ? width - 1 : width);

	for (int i = 0; i < iterations; ++i)
	{
		cout << "Descent iteration " << i + 1 << "..." << endl;

		for (int y = yStart; y < yEnd; ++y)
		{
			// The left and right boundaries are left untouched
			for (int x = xStart; x < xEnd; ++x)
			{
				auto xm = max(x - 1, 0);
				auto xp = min(x + 1, width - 1);
				auto ym = max(y - 1, 0);
				auto yp = min(y + 1, height - 1);
				auto df = laplacian[y][x] - ((*readField)[y][xp] + (*readField)[y][xm] + (*readField)[yp][x] + (*readField)[ym][x]
					- 4.0f * (*readField)[y][x]);
				(*writeField)[y][x] = fmax(fmin((*readField)[y][x] - epsilon * df, 1.0f), 0.0f);
			}
		}

		// Swap the fields
		auto temp = writeField;
		writeField = readField;
		readField = temp;
	}

	copyField(*readField, output);
}

// Solve using poisson equation
void poissonSolver(const ScalarField& margin, const ScalarField& boundary, ScalarField* output)
{
	assert(margin.size() == boundary.size());
	assert(margin[0].size() == boundary[0].size());
	resizeField(margin[0].size(), margin.size(), output);

	auto gridWidth = margin[0].size();
	auto gridHeight = margin.size();
	auto N = gridWidth * gridHeight;
	Eigen::MatrixXf A(N, N);
	Eigen::VectorXf b(N);
	Eigen::VectorXf z;
	unordered_map<int, pair<int, int>> index_to_gridindex;
	unordered_map<pair<int, int>, int, pairhash> gridindex_to_index;
	vector<bool> index_is_boundary;

	cout << "Setting up linear system..." << endl;

	for (size_t i = 0; i < N; ++i)
	{
		pair<int, int> gridindex;
		gridindex.first = i / gridWidth;
		gridindex.second = i % gridWidth;
		bool is_boundary = false;
		bool on_y_boundary = (gridindex.first == 0) || (gridindex.first == gridHeight - 1);
		bool on_x_boundary = (gridindex.second == 0) || (gridindex.second == gridWidth - 1);
		if (on_y_boundary || on_x_boundary)
		{
			is_boundary = true;
		}
		index_to_gridindex[i] = gridindex;
		gridindex_to_index[gridindex] = i;
		index_is_boundary.push_back(is_boundary);
	}

	for (size_t x = 0; x < N; ++x)
	{
		pair<int, int> gridindex = index_to_gridindex[x];
		for (size_t y = 0; y < N; ++y)
			A(x, y) = 0;

		if (index_is_boundary[x])
		{
			A(x, x) = 1;
			b(x) = boundary[gridindex.first][gridindex.second];
			continue;
		}

		b(x) = margin[gridindex.first][gridindex.second];
		unordered_set<int> active_variables;

		pair<int, int> gridindex_down(gridindex.first + 1, gridindex.second);
		pair<int, int> gridindex_up(gridindex.first - 1, gridindex.second);
		pair<int, int> gridindex_left(gridindex.first, gridindex.second - 1);
		pair<int, int> gridindex_right(gridindex.first, gridindex.second + 1);

		active_variables.insert(gridindex_to_index[gridindex_right]);
		active_variables.insert(gridindex_to_index[gridindex_left]);
		active_variables.insert(gridindex_to_index[gridindex_up]);
		active_variables.insert(gridindex_to_index[gridindex_down]);

		A(x, x) = -4;
		for (const auto& var : active_variables)
		{
			A(x, var) = 1;
		}
	}

	cout << "Solving system..." << endl;

	z = A.lu().solve(b);

	for (size_t i = 0; i < N; ++i)
	{
		pair<int, int> gridindex = index_to_gridindex[i];
		(*output)[gridindex.first][gridindex.second] = fmin(fmax(z(i), 0.0f), 1.0f);
	}
}

typedef Eigen::Triplet<double> DTriplet;
void insertCoefficient(int id, int x, int y, double w, const int width, const int height, std::vector<DTriplet>* coeffs,
	Eigen::VectorXd* b, const ScalarField& boundary)
{
	int targetVariableId = x + y * width;
	if (x == -1)
		(*b)(id) -= w * boundary[y + 1][0]; // constrained coefficient
	else if (x == width)
		(*b)(id) -= w * boundary[y + 1][width + 1]; // constrained coefficient
	else if (y == -1)
		(*b)(id) -= w * boundary[0][x + 1]; // constrained coefficient
	else if (y == height)
		(*b)(id) -= w * boundary[height + 1][x + 1]; // constrained coefficient
	else
		coeffs->push_back(DTriplet(id, targetVariableId, w)); // unknown coefficient
}

void buildProblem(std::vector<DTriplet>* coefficients, Eigen::VectorXd* b, const ScalarField& laplacian,
	const ScalarField& boundary)
{
	b->setZero();
	const int width = laplacian[0].size() - 2;
	const int height = laplacian.size() - 2;

	int id = 0;
	for (int y = 0; y < height; ++y)
	{
		for (int x = 0; x < width; ++x, ++id)
		{
			(*b)(id) += laplacian[y + 1][x + 1];

			// Discrete laplacian star
			insertCoefficient(id, x - 1, y, 1, width, height, coefficients, b, boundary);
			insertCoefficient(id, x + 1, y, 1, width, height, coefficients, b, boundary);
			insertCoefficient(id, x, y - 1, 1, width, height, coefficients, b, boundary);
			insertCoefficient(id, x, y + 1, 1, width, height, coefficients, b, boundary);
			insertCoefficient(id, x, y, -4, width, height, coefficients, b, boundary);
		}
	}
}

void poissonSolverSparse(const ScalarField& margin, const ScalarField& boundary, ScalarField* output)
{
	assert(margin.size() == boundary.size());
	assert(margin[0].size() == boundary[0].size());

	const int width = margin[0].size();
	const int height = margin.size();
	const int variableCount = (width - 2) * (height - 2);

	resizeField(width, height, output);

	// Construct the linear system
	Eigen::VectorXd b(variableCount);
	Eigen::SparseMatrix<double> A(variableCount, variableCount);
	std::vector<DTriplet> coefficients;

	cout << "Building sparse linear system..." << endl;

	buildProblem(&coefficients, &b, margin, boundary);
	A.setFromTriplets(coefficients.begin(), coefficients.end());

	cout << "Solving linear system (simplicial cholesky)..." << endl;

	// Solve the system
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > chol(A);
	Eigen::VectorXd result = chol.solve(b);

	copyField(boundary, output);

	// Copy over output
	int id = 0;
	for (int y = 1; y < height - 1; ++y)
		for (int x = 1; x < width - 1; ++x, ++id)
			(*output)[y][x] = static_cast<float>(fmin(fmax(result(id), 0.0f), 1.0f));
}

int main(int argc, char** argv)
{
	// Read command line inputs if specified
	string imageSource1 = DEFAULT_IMAGE_SOURCE_1;
	string imageSource2 = DEFAULT_IMAGE_SOURCE_2;
	string outputPath = DEFAULT_IMAGE_OUTPUT;
	int stitchMargin = DEFAULT_STITCH_MARGIN;
	ProgramMode mode = DEFAULT_MODE;
	bool bUseFileOutput = false;

	// Read input parameters if needed
	if (argc >= 3)
	{
		imageSource1 = argv[1];
		imageSource2 = argv[2];
	}
	if (argc >= 4)
		stitchMargin = stoi(argv[3]);
	if (argc >= 5)
		mode = static_cast<ProgramMode>(stoi(argv[4]));
	if (argc >= 6)
	{
		cout << "File output specified!" << endl;
		outputPath = argv[5];
		bUseFileOutput = true;
	}

	// Open the PNG files for image 1 and image 2
	ScalarField imageArray1;
	ScalarField imageArray2;
	unsigned int error1 = floatMatrixFromPNG(imageSource1, &imageArray1);
	unsigned int error2 = floatMatrixFromPNG(imageSource2, &imageArray2);

	if (error1 || error2)
	{
		cout << "Failed to open image files!" << endl;
		return EXIT_FAILURE;
	}

	// Resize images if necessary
	if (imageArray1.size() != imageArray2.size())
	{
		cout << "Unequal image height... resizing..." << endl;
		auto minsize = min(imageArray1.size(), imageArray2.size());
		if (imageArray1.size() != minsize)
			imageArray1.resize(minsize);
		else
			imageArray2.resize(minsize);
	}

	// Check margin
	if (stitchMargin < 1)
	{
		cout << "Invalid margin size!" << endl;
		return EXIT_FAILURE;
	}

	ScalarField output;
	bool bSkipSave = false;

	switch (mode)
	{
	case SimpleStitch:
	{
		// Stitch the images together
		ScalarField image1Margin;
		ScalarField image1Remainder;
		ScalarField image2Margin;
		ScalarField image2Remainder;
		subRange(imageArray1, 0, 0, imageArray1[0].size() - stitchMargin, imageArray1.size(), &image1Remainder);
		subRange(imageArray1, imageArray1[0].size() - stitchMargin, 0, stitchMargin, imageArray1.size(), &image1Margin);
		subRange(imageArray2, 0, 0, stitchMargin, imageArray2.size(), &image2Margin);
		subRange(imageArray2, stitchMargin, 0, imageArray2[0].size() - stitchMargin, imageArray2.size(), &image2Remainder);

		cout << "Stitching images..." << endl;
		ScalarField maskResult;
		computeStitchMask(image1Margin, image2Margin, &maskResult);
		ScalarField mergeResult;
		applyMask(image1Margin, image2Margin, maskResult, &mergeResult);
		cout << "Stitching complete!" << endl;

		cout << "Merging..." << endl;
		vector<ScalarField*> mergeParams = { &image1Remainder, &mergeResult, &image2Remainder };
		mergeRows(mergeParams, &output);
		cout << "Merging Complete!" << endl;

		break;
	}
	case RecoverFromGradient:
	{
		// Atempt to recover an image from gradient information
		VectorField gradient;
		computeGradient(imageArray1, &gradient);

		int width = imageArray1[0].size();
		int height = imageArray1.size();
		ScalarField field;
		resizeField(width, height, &field);

		for (int y = 0; y < height; ++y)
			for (int x = 0; x < width; ++x)
				field[y][x] = 0.5f;

		gradientDescent(field, gradient, RECOVERY_GRADIENT_ITERATIONS,
			GRADIENT_DESCENT_EPSILON, false, false, &output);

		break;
	}
	case RecoverFromLaplace:
	{
		ScalarField laplacian;
		computeLaplacian(imageArray1, &laplacian);

		int width = imageArray1[0].size();
		int height = imageArray1.size();
		ScalarField field;
		resizeField(width, height, &field);

		for (int y = 0; y < height; ++y)
			for (int x = 0; x < width; ++x)
				field[y][x] = 0.5f;

		laplacianDescent(field, laplacian, RECOVERY_LAPLACE_ITERATIONS,
			LAPLACIAN_DESCENT_EPSILON, false, false, &output);

		break;
	}
	case GradientStitch:
	{
		// Stitch the images together
		ScalarField image1Margin;
		ScalarField image1Remainder;
		ScalarField image2Margin;
		ScalarField image2Remainder;
		subRange(imageArray1, 0, 0, imageArray1[0].size() - stitchMargin, imageArray1.size(), &image1Remainder);
		subRange(imageArray1, imageArray1[0].size() - stitchMargin, 0, stitchMargin, imageArray1.size(), &image1Margin);
		subRange(imageArray2, 0, 0, stitchMargin, imageArray2.size(), &image2Margin);
		subRange(imageArray2, stitchMargin, 0, imageArray2[0].size() - stitchMargin, imageArray2.size(), &image2Remainder);

		// Stitch the gradients together
		VectorField marginGradient1;
		VectorField marginGradient2;
		computeGradient(image1Margin, &marginGradient1);
		computeGradient(image2Margin, &marginGradient2);

		ScalarField mask;
		cout << "Computing mask..." << endl;
		computeGradientStitchMask(marginGradient1, marginGradient2, &mask);
		cout << "Mask computation complete!" << endl;

		// Perform gradient descent
		VectorField mergedGradient;
		applyMask(marginGradient1, marginGradient2, mask, &mergedGradient);
		ScalarField initialGuess;
		applyMask(image1Margin, image2Margin, mask, &initialGuess);
		ScalarField result;
		gradientDescent(initialGuess, mergedGradient, GRADIENT_DESCENT_ITERATIONS,
			GRADIENT_DESCENT_EPSILON, true, false, &result);

		cout << "Merging..." << endl;
		vector<ScalarField*> mergeParams = { &image1Remainder, &result, &image2Remainder };
		mergeRows(mergeParams, &output);
		cout << "Merging Complete!" << endl;
		break;
	}
	case LaplaceStitch:
	{
		// Stitch the images together
		ScalarField image1Margin;
		ScalarField image1Remainder;
		ScalarField image2Margin;
		ScalarField image2Remainder;
		subRange(imageArray1, 0, 0, imageArray1[0].size() - stitchMargin, imageArray1.size(), &image1Remainder);
		subRange(imageArray1, imageArray1[0].size() - stitchMargin, 0, stitchMargin, imageArray1.size(), &image1Margin);
		subRange(imageArray2, 0, 0, stitchMargin, imageArray2.size(), &image2Margin);
		subRange(imageArray2, stitchMargin, 0, imageArray2[0].size() - stitchMargin, imageArray2.size(), &image2Remainder);

		// Stitch the gradients together
		VectorField marginGradient1;
		VectorField marginGradient2;
		computeGradient(image1Margin, &marginGradient1);
		computeGradient(image2Margin, &marginGradient2);

		ScalarField mask;
		cout << "Computing mask..." << endl;
		computeGradientStitchMask(marginGradient1, marginGradient2, &mask);
		cout << "Mask computation complete!" << endl;

		cout << "Computing laplacian..." << endl;
		ScalarField marginLaplacian1;
		ScalarField marginLaplacian2;
		computeLaplacian(image1Margin, &marginLaplacian1);
		computeLaplacian(image2Margin, &marginLaplacian2);
		ScalarField laplacian;
		ScalarField initialGuess;
		applyMask(marginLaplacian1, marginLaplacian2, mask, &laplacian);
		applyMask(image1Margin, image2Margin, mask, &initialGuess);

		// Perform laplacian descent
		ScalarField result;
		laplacianDescent(initialGuess, laplacian, LAPLACIAN_DESCENT_ITERATIONS,
			LAPLACIAN_DESCENT_EPSILON, true, false, &result);

		cout << "Merging..." << endl;
		vector<ScalarField*> mergeParams = { &image1Remainder, &result, &image2Remainder };
		mergeRows(mergeParams, &output);
		cout << "Merging Complete!" << endl;
		break;
	}
	case PoissonStitch:
	case PoissonStitchSparse:
	{
		if (mode == PoissonStitch)
			cout << "Warning: Sparse poisson mode is recommended for large images..." << endl;

		// Stitch the images together
		ScalarField image1Margin;
		ScalarField image1Remainder;
		ScalarField image2Margin;
		ScalarField image2Remainder;
		subRange(imageArray1, 0, 0, imageArray1[0].size() - stitchMargin, imageArray1.size(), &image1Remainder);
		subRange(imageArray1, imageArray1[0].size() - stitchMargin, 0, stitchMargin, imageArray1.size(), &image1Margin);
		subRange(imageArray2, 0, 0, stitchMargin, imageArray2.size(), &image2Margin);
		subRange(imageArray2, stitchMargin, 0, imageArray2[0].size() - stitchMargin, imageArray2.size(), &image2Remainder);

		// Stitch the gradients together
		VectorField marginGradient1;
		VectorField marginGradient2;
		computeGradient(image1Margin, &marginGradient1);
		computeGradient(image2Margin, &marginGradient2);

		ScalarField mask;
		cout << "Computing mask..." << endl;
		computeGradientStitchMask(marginGradient1, marginGradient2, &mask);
		cout << "Mask computation complete!" << endl;

		cout << "Computing laplacian..." << endl;
		ScalarField marginLaplacian1;
		ScalarField marginLaplacian2;
		computeLaplacian(image1Margin, &marginLaplacian1);
		computeLaplacian(image2Margin, &marginLaplacian2);
		ScalarField laplacian;
		ScalarField boundary;
		applyMask(marginLaplacian1, marginLaplacian2, mask, &laplacian);
		applyMask(image1Margin, image2Margin, mask, &boundary);

		cout << "Running Poisson solver..." << endl;
		ScalarField result;
		if (mode == PoissonStitch)
			poissonSolver(laplacian, boundary, &result);
		else if (mode == PoissonStitchSparse)
			poissonSolverSparse(laplacian, boundary, &result);

		cout << "Merging..." << endl;
		vector<ScalarField*> mergeParams = { &image1Remainder, &result, &image2Remainder };
		mergeRows(mergeParams, &output);
		cout << "Merging Complete!" << endl;
		break;
	}
	case ComputeGradient:
	{
		// Compute the gradient of an image and save the result
		VectorField gradient;
		computeGradient(imageArray1, &gradient);
		vector<unsigned char> outputData;
		convertGradientToImageData(gradient, &outputData);
		auto error = lodepng::encode(outputPath, outputData, static_cast<unsigned int>(gradient[0].size()),
			static_cast<unsigned int>(gradient.size()));

		bSkipSave = true;
		break;
	}
	default:
		cout << "Unrecognized command!" << endl;
		return EXIT_SUCCESS;
	}

	if (!bSkipSave)
	{
		// Display mode
		if (!bUseFileOutput)
		{
#ifdef _WIN32
			// Display the output in an opengl window (windows-only)
			cout << "Displaying result..." << endl;
			vector<unsigned char> imageData;
			convertFloatMatrixToImageData(output, &imageData);

			displayImage(output[0].size(), output.size(), imageData.data());
#else
			// No display mode on non-windows platforms
			bUseFileOutput = true;
#endif
		}
		
		// File output mode
		if (bUseFileOutput)
		{
			// Save the result
			cout << "Saving result..." << endl;
			auto error = saveFloatMatrixToPNG(outputPath, output);

			if (error)
			{
				cout << "Failed to save result!" << endl;
				return EXIT_FAILURE;
			}
		}
	}

	cout << "Success!" << endl;
	return EXIT_SUCCESS;
}