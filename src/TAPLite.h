#ifndef GUARD_PATH_ENGINE_H
#define GUARD_PATH_ENGINE_H

#ifdef _WIN32
#define PATH_ENGINE_API __declspec(dllexport)
#else
#define PATH_ENGINE_API
#endif

extern "C" PATH_ENGINE_API void DTA_AssignmentAPI();

extern "C" PATH_ENGINE_API void DTA_SimulationAPI();

#define BUFFERSIZE 1000
#define MAX_NO_BISECTITERATION 5 /* Avoids infinite loops */

#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

static int MinLineSearchIterations = 5;
static int ActualIterations = 0;
static double LastLambda = 1.0;
int g_accessibility_only_mode = 0;

std::map<int, int> g_map_external_node_id_2_node_seq_no;
std::map<int, int> g_map_node_seq_no_2_external_node_id;
std::map<int, int> g_map_internal_zone_no_2_node_seq_no;
class CDTACSVParser {
public:
    char Delimiter;
    bool IsFirstLineHeader;
    // for DataHub CSV files
    bool m_bSkipFirstLine;
    bool m_bDataHubSingleCSVFile;
    bool m_bLastSectionRead;

    std::ifstream inFile;
    std::string mFileName;
    std::string m_DataHubSectionName;
    std::string SectionName;

    std::vector<std::string> LineFieldsValue;
    std::vector<int> LineIntegerVector;
    std::vector<std::string> Headers;
    std::map<std::string, int> FieldsIndices;

    CDTACSVParser()
        : Delimiter{ ',' },
        IsFirstLineHeader{ true },
        m_bSkipFirstLine{ false },
        m_bDataHubSingleCSVFile{ false },
        m_bLastSectionRead{ false }
    {
    }

    ~CDTACSVParser()
    {
        if (inFile.is_open())
            inFile.close();
    }

    // inline member functions
    std::vector<std::string> GetHeaderVector() { return Headers; }
    void CloseCSVFile() { inFile.close(); }

    void ConvertLineStringValueToIntegers();
    bool OpenCSVFile(std::string fileName, bool b_required);
    bool ReadRecord();
    bool ReadSectionHeader(std::string s);
    bool ReadRecord_Section();
    std::vector<std::string> ParseLine(std::string line);
    bool GetValueByFieldName(std::string field_name,
        std::string& value,
        bool required_field = true);
    template <class T>
    bool GetValueByFieldName(std::string field_name,
        T& value,
        bool required_field = true,
        bool NonnegativeFlag = true);
    template <class T>
    bool GetValueByKeyName(std::string field_name,
        T& value,
        bool required_field = true,
        bool NonnegativeFlag = true);
};

// definitions of CDTACSVParser member functions
void CDTACSVParser::ConvertLineStringValueToIntegers()
{
    LineIntegerVector.clear();
    for (unsigned i = 0; i < LineFieldsValue.size(); ++i)
    {
        std::string si = LineFieldsValue[i];
        int value = atoi(si.c_str());

        if (value >= 1)
            LineIntegerVector.push_back(value);
    }
}

// Trim in place
static inline void TrimString(std::string& s)
{
    // Left trim
    s.erase(
        s.begin(),
        std::find_if(s.begin(), s.end(), [](unsigned char ch) { return !std::isspace(ch); })
    );

    // Right trim
    s.erase(
        std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(),
        s.end()
    );
}

bool CDTACSVParser::OpenCSVFile(std::string fileName, bool b_required)
{
    mFileName = fileName;
    inFile.open(fileName.c_str());

    if (!inFile.is_open())
    {
        // If this CSV file is required and cannot be opened, handle error
        if (b_required)
        {
            // g_program_stop(); // or throw an exception, or however you prefer
        }
        return false;
    }

    // If the first line is a header row, read and parse it
    if (IsFirstLineHeader)
    {
        std::string s;
        if (!std::getline(inFile, s))
        {
            // File is open but no header line found
            return false;
        }

        // -- Optional: Remove UTF-8 BOM if present --
        static const unsigned char UTF8_BOM[3] = { 0xEF, 0xBB, 0xBF };
        if (s.size() >= 3 &&
            (unsigned char)s[0] == UTF8_BOM[0] &&
            (unsigned char)s[1] == UTF8_BOM[1] &&
            (unsigned char)s[2] == UTF8_BOM[2])
        {
            s.erase(0, 3);
        }

        // Parse the first line to get header fields
        std::vector<std::string> FieldNames = ParseLine(s);

        // Trim each header field and store in FieldsIndices map
        for (size_t i = 0; i < FieldNames.size(); ++i)
        {
            std::string tmp_str = FieldNames[i];
            TrimString(tmp_str);

            // You might decide to do further normalization here, e.g. transform to lowercase
            // std::transform(tmp_str.begin(), tmp_str.end(), tmp_str.begin(), ::tolower);

            FieldsIndices[tmp_str] = static_cast<int>(i);
        }
    }

    return true;
}
bool CDTACSVParser::ReadRecord()
{
    LineFieldsValue.clear();

    if (inFile.is_open())
    {
        std::string s;
        std::getline(inFile, s);
        if (s.length() > 0)
        {
            LineFieldsValue = ParseLine(s);
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

std::vector<std::string> CDTACSVParser::ParseLine(std::string line)
{
    std::vector<std::string> SeperatedStrings;
    std::string subStr;

    if (line.length() == 0)
        return SeperatedStrings;

    std::istringstream ss(line);

    if (line.find_first_of('"') == std::string::npos)
    {
        while (std::getline(ss, subStr, Delimiter))
        {
            SeperatedStrings.push_back(subStr);
        }

        if (line.at(line.length() - 1) == ',')
        {
            SeperatedStrings.push_back("");
        }
    }
    else
    {
        while (line.length() > 0)
        {
            size_t n1 = line.find_first_of(',');
            size_t n2 = line.find_first_of('"');

            if (n1 == std::string::npos &&
                n2 == std::string::npos)  // last field without double quotes
            {
                subStr = line;
                SeperatedStrings.push_back(subStr);
                break;
            }

            if (n1 == std::string::npos && n2 != std::string::npos)  // last field with double
                // quotes
            {
                size_t n3 = line.find_first_of('"', n2 + 1);  // second double quote

                // extract content from double quotes
                subStr = line.substr(n2 + 1, n3 - n2 - 1);
                SeperatedStrings.push_back(subStr);

                break;
            }

            if (n1 != std::string::npos && (n1 < n2 || n2 == std::string::npos))
            {
                subStr = line.substr(0, n1);
                SeperatedStrings.push_back(subStr);
                if (n1 < line.length() - 1)
                {
                    line = line.substr(n1 + 1);
                }
                else  // comma is the last char in the line string, push an empty string to the back
                      // of vector
                {
                    SeperatedStrings.push_back("");
                    break;
                }
            }

            if (n1 != std::string::npos && n2 != std::string::npos && n2 < n1)
            {
                size_t n3 = line.find_first_of('"', n2 + 1);  // second double quote
                subStr = line.substr(n2 + 1, n3 - n2 - 1);
                SeperatedStrings.push_back(subStr);
                size_t idx = line.find_first_of(',', n3 + 1);

                if (idx != std::string::npos)
                {
                    line = line.substr(idx + 1);
                }
                else
                {
                    break;
                }
            }
        }
    }
    return SeperatedStrings;
}

bool CDTACSVParser::GetValueByFieldName(std::string field_name,
    std::string& value,
    bool required_field)
{
    if (FieldsIndices.find(field_name) == FieldsIndices.end())
    {
        if (required_field)
        {
            // dtalog.output() << "[ERROR] Field " << field_name << " in file " << mFileName << "
            // does not exist. Please check the file." << '\n'; g_DTA_log_file << "[ERROR] Field "
            // << field_name << " in file " << mFileName << " does not exist. Please check the
            // file." << '\n'; g_program_stop();
        }
        return false;
    }
    else
    {
        if (LineFieldsValue.size() == 0)
        {
            return false;
        }

        unsigned int index = FieldsIndices[field_name];
        if (index >= LineFieldsValue.size())
        {
            return false;
        }
        std::string str_value = LineFieldsValue[index];

        if (str_value.length() <= 0)
        {
            return false;
        }

        value = str_value;
        return true;
    }
}

// Peiheng, 03/22/21, to avoid implicit instantiations in flash_dta.cpp and main_api.cpp for this
// template function only all the other non-inline functions are implemented in utils.cpp
template <class T>
bool CDTACSVParser::GetValueByFieldName(std::string field_name,
    T& value,
    bool required_field,
    bool NonnegativeFlag)
{
    if (FieldsIndices.find(field_name) == FieldsIndices.end())
    {
        if (required_field)
        {
            //cout << "[ERROR] Field " << field_name << " in file " << mFileName.c_str()
            // << " does not exist. Please check the file." << '\n'; 
        }
        return false;
    }
    else
    {
        if (LineFieldsValue.size() == 0)
        {
            return false;
        }

        int size = (int)(LineFieldsValue.size());
        if (FieldsIndices[field_name] >= size)
        {
            return false;
        }

        std::string str_value = LineFieldsValue[FieldsIndices[field_name]];

        if (str_value.length() <= 0)
        {
            return false;
        }

        std::istringstream ss(str_value);

        T converted_value;
        ss >> converted_value;

        if (/*!ss.eof() || */ ss.fail())
        {
            return false;
        }

        // if (required_field)
        //{
        //     if(NonnegativeFlag)
        //     {
        //         if (converted_value < 0)
        //             converted_value = 0;
        //     }
        // }

        value = converted_value;
        return true;
    }
}

void ExitMessage(const char* format, ...)
{
	va_list ap;

	// vprintf(format, ap);
	printf("\n");

	getchar();

	exit(EXIT_FAILURE);
}

// Function to allocate a 1D array
void* Alloc_1D(int dim1, size_t size) {
    void* Array = (void*)calloc(dim1 + 1, size);
    if (Array == NULL) {
        ExitMessage("Cannot allocate memory for single dimension array of size %d.\n", dim1);
    }
    return Array;
}

// Function to allocate a 2D array
void** Alloc_2D(int dim1, int dim2, size_t size) {
    void** Array = (void**)calloc(dim1 + 1, sizeof(void*));
    int i;  // Loop variable declared outside for efficiency

    if (Array == NULL) {
        ExitMessage("Cannot allocate memory for two-dimensional array of size %d by %d.\n", dim1, dim2);
    }

    for (i = 0; i <= dim1; i++) {
        Array[i] = (void*)calloc(dim2 + 1, size);
        if (Array[i] == NULL) {
            ExitMessage("Cannot allocate memory for two-dimensional array of size %d by %d.\n", dim1, dim2);
        }
    }
    return Array;
}

// Function to free a 2D array
void Free_2D(void** Array, int dim1, int dim2) {
    int i;  // Loop variable declared outside for efficiency

    for (i = 0; i <= dim1; i++) {
        free(Array[i]);
    }
    free(Array);
}

// Function to allocate a 3D array
void*** Alloc_3D(int dim1, int dim2, int dim3, size_t size) {
    void*** Array = (void***)calloc(dim1 + 1, sizeof(void**));
    int i, j;  // Loop variables declared outside for efficiency

    if (Array == NULL) {
        ExitMessage("Cannot allocate memory for three-dimensional array of size %d by %d by %d.\n", dim1, dim2, dim3);
    }

    for (i = 0; i <= dim1; i++) {
        Array[i] = (void**)calloc(dim2 + 1, sizeof(void*));
        if (Array[i] == NULL) {
            ExitMessage("Cannot allocate memory for two-dimensional array of size %d by %d.\n", dim1, dim2);
        }

        for (j = 0; j <= dim2; j++) {
            Array[i][j] = (void*)calloc(dim3 + 1, size);
            if (Array[i][j] == NULL) {
                ExitMessage("Cannot allocate memory for one-dimensional array of size %d.\n", dim3);
            }
        }
    }
    return Array;
}

// Function to free a 3D array
void Free_3D(void*** Array, int dim1, int dim2, int dim3) {
    int i, j;   // Loop variables declared outside for efficiency
    void* p;    // Pointer variable declared outside the loop for efficiency

    // Free the innermost arrays (1D arrays)
    for (i = 0; i <= dim1; i++) {
        for (j = 0; j < dim2; j++) {
            p = Array[i][j];
            free(p);
        }
        // Free the 2D array (array of pointers to 1D arrays)
        free(Array[i]);
    }

    // Free the outermost array (array of pointers to 2D arrays)
    free(Array);
}

/* Internal functions */

/* Initialization functions */
struct CLink {
    int link_id;
    int internal_from_node_id;
    int internal_to_node_id;
    int length;
    int lanes;
    double capacity;
    int free_speed;
};

#endif
