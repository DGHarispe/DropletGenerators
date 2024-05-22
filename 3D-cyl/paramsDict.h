#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_DICTIONARY_SIZE 100

struct KeyValuePair { const char* key; const char* value;};
struct KeyValuePair dictionary[100] = {
    {"NCells", "int"},
	{"maxRlevel", "int"},
	{"q1", "float"},
    {"q2", "float"},
	{"SIGMA", "float"},
	{"radius", "float"},
	{"uemax", "float"},
	{"ownL0", "float"},
	{"size_domm", "float"},
	{"rho1_val", "float"},
	{"rho2_val", "float"},
	{"mu1_val", "float"},
	{"mu2_val", "float"},
	{"oxc", "float"},
	{"oyc", "float"},
	{"ozc", "float"},
	{"DTi", "float"},
	{"nMaxDT", "float"},
	{"startSteps", "int"},
	{"tmax", "float"},
	{"csTol", "float"},
	{"fTol", "float"}
};

int   	NCells, maxRlevel;
float	q1, q2;
float	SIGMA, radius, uemax, ownL0, size_domm, rho1_val, rho2_val, mu1_val, mu2_val;
float	oxc, oyc, ozc;
float	DTi, nMaxDT;
int		startSteps;
float	tmax, csTol, fTol;

void* variableArray[] = 
{
	&NCells, &maxRlevel,
	&q1, &q2, &SIGMA, &radius, &uemax, &ownL0, &size_domm, &rho1_val, &rho2_val, &mu1_val, &mu2_val,
	&oxc, &oyc, &ozc, &DTi, &nMaxDT, &startSteps, &tmax, &csTol, &fTol
};

const char* getValue(const char* key) {
    for (int i = 0; i < MAX_DICTIONARY_SIZE; i++) {
        if (strcmp(key, dictionary[i].key) == 0) {
            return dictionary[i].value;
        }
    }
    return NULL;
}

void setValue(const char* key, const char* value) {
    for (int i = 0; i < MAX_DICTIONARY_SIZE; i++) {
        if (dictionary[i].key == NULL) {
            dictionary[i].key = key;
            dictionary[i].value = value;
            return;
        }
    }
    printf("Dictionary full!\n");
}

void setParamValues(void* variable, char* key, char *value) {
	const char* varType = getValue(key);
	
    if (strcmp(varType, "int") == 0) {
		*(int*) variable = atoi(value);
    }
    else if (strcmp(varType, "float") == 0) {
		*(float*) variable = atof(value);
    }
    else if (strcmp(varType, "string") == 0) {
		int length = strlen(value);
		strncpy(variable, value, length);
		//variable[sizeof(variable) - 1] = '\0';  // Ensure null-termination
    }
    else {
        printf("Invalid param type!: %s\n", getValue(key));
    }
}

bool readParams() {
    FILE *file;
    char line[200];
    char *param, *value;

    file = fopen("parameters.txt", "r");
    if (file == NULL) {
        printf("Failed to open the file.\n");
        return false;
    }

    for (int varIndex = 0; fgets(line, sizeof(line), file); varIndex++) {
        if (line[0] == '#'){	//Ignore comments
			varIndex--;
            continue; }
		
        param = strtok(line, ": ");
        value = strtok(NULL, "");
		
        if (param == NULL || value == NULL){
			varIndex--;
			continue; }
		
		// Remove leading/trailing whitespaces from the value
		value[strcspn(value, "\r\n")] = '\0';
		
		setParamValues(variableArray[varIndex], param, value);
    }
	
    fclose(file);
    return true;
}