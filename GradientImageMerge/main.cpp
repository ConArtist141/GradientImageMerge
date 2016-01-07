#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>
#include <string>
#include <algorithm>
#include <cmath>

#include "lodepng.h"
#include "GridGraph_2D_4C.h"

#define DEFAULT_IMAGE_SOURCE_1 "cat.png"
#define DEFAULT_IMAGE_SOURCE_2 "cat.png"
#define DEFAULT_IMAGE_OUTPUT "result.png"
#define DEFAULT_STITCH_MARGIN 100
#define DEFAULT_MODE SimpleStitch

using namespace std;

// Program mode
enum ProgramMode
{
	SimpleStitch,
	ComputeGradient,
	GradientStitch
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

typedef vector<vector<float> > ScalarField;
typedef vector<vector<vec2<float> > > VectorField;
typedef GridGraph_2D_4C<float, float, float> Grid;

void performStitching(ScalarField field1, ScalarField field2, ScalarField* mergeResult)
{
	assert(field1.size() == field2.size());
	assert(field1[0].size() == field2[0].size());

	auto gridWidth = field1[0].size();
	auto gridHeight = field1.size();

	mergeResult->resize(gridHeight);
	for (auto& row : *mergeResult)
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
				(*mergeResult)[y][x] = field1[y][x];
			else
				(*mergeResult)[y][x] = field2[y][x];
		}
	}
}

template <typename t>
void mergeRows(const vector<vector<vector<t> >* >& fields, vector<vector<t> >* output)
{
	auto fieldWidth = 0;
	auto fieldHeight = (*fields[0])[0].size();
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

// Convert gradient data to image data
void convertGradientToImageData(VectorField& grad, vector<unsigned char>* output)
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

int main(int argc, char** argv)
{
	// Read command line inputs if specified
	string imageSource1 = DEFAULT_IMAGE_SOURCE_1;
	string imageSource2 = DEFAULT_IMAGE_SOURCE_2;
	string outputPath = DEFAULT_IMAGE_OUTPUT;
	int stitchMargin = DEFAULT_STITCH_MARGIN;
	ProgramMode mode = DEFAULT_MODE;

	// Read input parameters if needed
	if (argc >= 2)
	{
		imageSource1 = argv[0];
		imageSource2 = argv[1];
	}
	if (argc >= 3)
		stitchMargin = stoi(argv[2]);
	if (argc >= 4)
		outputPath = argv[3];
	if (argc >= 5)
		mode = static_cast<ProgramMode>(stoi(argv[4]));

	// Open the PNG files for image 1 and image 2
	ScalarField imageArray1;
	ScalarField imageArray2;
	unsigned int error1 = floatMatrixFromPNG(imageSource1, &imageArray1);
	unsigned int error2 = floatMatrixFromPNG(imageSource2, &imageArray2);

	if (error1 || error2)
	{
		cout << "Failed to open image files!" << endl;
		return -1;
	}

	unsigned int error;
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
		ScalarField mergeResult;
		performStitching(image1Margin, image2Margin, &mergeResult);
		cout << "Stitching complete!" << endl;

		cout << "Merging..." << endl;
		vector<ScalarField*> mergeParams = { &image1Remainder, &mergeResult, &image2Remainder };
		ScalarField output;
		mergeRows(mergeParams, &output);
		cout << "Merging Complete!" << endl;

		// Save the result
		cout << "Saving result..." << endl;
		error = saveFloatMatrixToPNG(outputPath, output);
		break;
	}
	case ComputeGradient:
	{
		// Compute the gradient of an image and save the result
		VectorField gradient;
		computeGradient(imageArray1, &gradient);
		vector<unsigned char> outputData;
		convertGradientToImageData(gradient, &outputData);
		error = lodepng::encode(outputPath, outputData, static_cast<unsigned int>(gradient[0].size()),
			static_cast<unsigned int>(gradient.size()));
		break;
	}
	case GradientStitch:
	{
		// Stitch the gradients together
		VectorField gradient1;
		VectorField gradient2;
		computeGradient(imageArray1, &gradient1);
		computeGradient(imageArray2, &gradient2);

		VectorField gradientOutput;
		cout << "Stitching gradients..." << endl;
		// performGradientStitching(gradient1, gradient2, stitchMargin, &gradientOutput);
		cout << "Stitching complete!" << endl;

		// Save the result
		cout << "Saving result..." << endl;
		vector<unsigned char> outputData;
		convertGradientToImageData(gradientOutput, &outputData);
		error = lodepng::encode(outputPath, outputData, static_cast<unsigned int>(gradientOutput[0].size()),
			static_cast<unsigned int>(gradientOutput.size()));
		break;
	}
	}

	if (error)
	{
		cout << "Failed to save result!" << endl;
		return -1;
	}

	cout << "Success!" << endl;
	return 0;
}