#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <stdio.h>
#include <stdlib.h>  
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <queue>
#define DEBUG 0
using namespace std;

typedef struct Pin {
    string name;
    string direct;
    double cap;
} Pin;

typedef struct CellLibrary {
    string name;
    map<string, Pin> pins;
    vector<double> index_1, index_2;
    vector< vector<double> > cell_rise;
    vector< vector<double> > cell_fall;
    vector< vector<double> > rise_transition;
    vector< vector<double> > fall_transition;
} CellLibrary;

typedef struct Gate {
    string name;
    string type;
    vector<string> ports;
    vector<string> wires;
    vector< vector<bool> > wires_input;
    vector<bool> output;
    int input_num_in_topological;
    bool isOutput;
    list<int> next;
    list<int> prev;
    bool input;
    vector<double> critical_cell_delay;
    vector<double> critical_transition;
} Gate;

class Module {
  public:
    string name;
    vector<string> input;
    vector<string> output;
    vector<string> wire;
    vector< vector<bool> > input_signal;
    vector< vector<bool> > output_signal;
    vector<Gate> gates;
    vector<int> output_gates;
    vector<int> topological_order;
    vector< vector<double> > max_delay;
    vector< vector<int> > max_delay_path_from_prev_gate;
    int critical_path_out_index;
    void topological_sort();
    void delay(map<string, CellLibrary> &cells);
    void output_file();

  private:
    void output_caculate();
    void max_delay_path_calculation(int index, double total_output_net_cap, CellLibrary &cell);
    double lookup_table(double &x, double &y, vector<double> &x_axis, vector<double> &y_axis, vector< vector<double> > &table);
};

void Module::output_caculate(){

    vector<bool> input_sig;
    vector<bool> output_sig;
    for(int x = 0; x < input_signal.size(); x++){ // patten number
        for(int i = 0; i < topological_order.size(); ++i){
            int idx = topological_order[i];
            for(int k = 0; k < gates[idx].wires.size()-1; k++){
                for(int j = 0; j < input.size(); j++){
                    if(gates[idx].wires[k] == input[j]){ // gate input == INPUT signal
                        input_sig.push_back(input_signal[x][j]);
                        break;
                    }
                }
            }
            if(input_sig.size()!=gates[idx].wires.size()-1){ // gate input  == prev gate output
                for (list<int>::iterator i = gates[idx].prev.begin(); i != gates[idx].prev.end(); ++i){ 
                    input_sig.push_back(gates[*i].output[x]);
                }
            }
            if(gates[idx].type == "NOR2X1"){
                if(!input_sig[0] && !input_sig[1]){
                    gates[idx].output.push_back(true);
                }
                else{
                    gates[idx].output.push_back(false);
                }
                #ifdef DEBUG
                    // cout << "----------------" << endl;
                    // cout << "gates : " << gates[idx].name << endl;
                    // cout << "Input : ";
                    // for(int i =0; i < input_sig.size(); i++){
                    //     cout << input_sig[i] << ", ";
                    // }
                    // cout << endl;
                    // cout << "Output : ";
                    // for(int i =0; i < gates[idx].output.size(); i++){
                    //     cout << gates[idx].output[i] << ", ";
                    // }
                    // cout << endl;
                    // cout << "----------------" << endl;
                #endif
                gates[idx].wires_input.push_back(input_sig);
                input_sig.clear();
            }
            else if(gates[idx].type == "INVX1"){
                if(input_sig[0]){
                    gates[idx].output.push_back(false);
                }
                else{
                    gates[idx].output.push_back(true);
                }
                #ifdef DEBUG
                // cout << "----------------" << endl;
                // cout << "gates : " << gates[idx].name << endl;
                // cout << "Input : ";
                // for(int i =0; i < input_sig.size(); i++){
                //     cout << input_sig[i] << ", ";
                // }
                // cout << endl;
                // cout << "Output : ";
                // for(int i =0; i < gates[idx].output.size(); i++){
                //     cout << gates[idx].output[i] << ", ";
                // }
                // cout << endl;
                // cout << "----------------" << endl;
                #endif
                gates[idx].wires_input.push_back(input_sig);
                input_sig.clear();
            }
            else if(gates[idx].type == "NANDX1"){
                if(input_sig[0] && input_sig[1]){
                    gates[idx].output.push_back(false);
                }
                else{
                    gates[idx].output.push_back(true);
                }
                #ifdef DEBUG
                // cout << "----------------" << endl;
                // cout << "gates : " << gates[idx].name << endl;
                // cout << "Input : ";
                // for(int i =0; i < input_sig.size(); i++){
                //     cout << input_sig[i] << ", ";
                // }
                // cout << endl;
                // cout << "Output : ";
                // for(int i =0; i < gates[idx].output.size(); i++){
                //     cout << gates[idx].output[i] << ", ";
                // }
                // cout << endl;
                // cout << "----------------" << endl;
                #endif
                gates[idx].wires_input.push_back(input_sig);
                input_sig.clear();
            }
        }
    }
}

void Module::topological_sort(){
	queue<int> q;
	for(int i = 0; i < gates.size(); i++){
		gates[i].input_num_in_topological = 0;
		gates[i].isOutput = false;
        // find output gate
		for(int j = 0; j < output.size(); j++){
            if(output[j] == gates[i].wires.back()){
                output_gates.push_back(i);
				gates[i].isOutput = true;
				break;
			}
		}
        // find input gate
		for(int j = 0; j < input.size(); j++){
			for(int k = 0; k < gates[i].wires.size()-1; k++){
                if(input[j] == gates[i].wires[k]){
					gates[i].input_num_in_topological++;
					break;
				}
			}
		}
		if(gates[i].input_num_in_topological == gates[i].wires.size()-1){ 
            q.push(i);
		}
	}
	while(!q.empty()){
		int idx;
        string first_gate_out;
        idx = q.front();
        q.pop();
        topological_order.push_back(idx);
		first_gate_out = gates[idx].wires.back();
		for(int i = 0; i < gates.size(); i++){
			for(int j = 0; j < gates[i].wires.size()-1; j++){
				if(gates[i].wires[j] == first_gate_out){ 
					gates[i].input_num_in_topological++;
					gates[i].prev.push_back(idx);
					gates[idx].next.push_back(i);
					if(gates[i].input_num_in_topological == gates[i].wires.size()-1){
                        q.push(i);
					}
					break;
				}
			}
		}
	}
    output_caculate();
// #ifdef DEBUG
//     cout << "--------------------------" << endl;
//     cout << "topological_order:" << endl;
//     for(int i = 0; i < topological_order.size(); ++i)
//         cout << topological_order[i] <<  ", ";
//     cout << endl;
//     cout << "--------------------------" << endl;
// #endif
#ifdef DEBUG
    cout << "--------------------------" << endl;
    cout << "topological result:" << endl;
    for(int i = 0; i < topological_order.size(); ++i){
        int idx = topological_order[i];
        cout << "gate No." << idx << "  " << gates[idx].name << endl;
        cout << "\tprev gates: ";
        for (list<int>::iterator i = gates[idx].prev.begin(); i != gates[idx].prev.end(); ++i){ 
            cout << *i << "(" << gates[*i].name << ")"<< ", ";
        }
        cout << endl;
        cout << "\tnext gates: ";
        for (list<int>::iterator j = gates[idx].next.begin(); j != gates[idx].next.end(); ++j){
            cout << *j << "(" << gates[*j].name << ")" << ", ";
        }
        cout << endl;
        cout << "\tInput boolean: " << endl;
        vector< vector<bool> >::iterator item2v = gates[idx].wires_input.begin();
        while(item2v !=gates[idx].wires_input.end()){
            vector<bool>::iterator item = (*item2v).begin();
            while(item != (*item2v).end()){
                cout<< "\t" << *item<<" , ";
                item++;
            }
            cout << endl;
            item2v++;
        }
        cout << "\tOutput boolean: ";
        for(int i =0; i < gates[idx].output.size(); i++){
            cout << gates[idx].output[i] << ", ";
        }
        cout << endl;
    }
    cout << "--------------------------" << endl;
#endif
}

double Module::lookup_table(double &row, double &col, vector<double> &row_index, vector<double> &col_index, vector< vector<double> > &table){
	int row_end_index = row_index.size()-1;
	int col_end_index = col_index.size()-1;
	int index_row1 = row_end_index;
	int index_row2 = 0;
	int index_col1 = col_end_index; 
	int index_col2 = 0;

	for(int i = row_end_index; i >= 0; i--) if(row_index[i] >  row) index_row1--;
	for(int i = 0; i <= row_end_index; i++) if(row_index[i] <= row) index_row2++;
	for(int i = col_end_index; i >= 0; i--) if(col_index[i] >  col) index_col1--;
	for(int i = 0; i <= col_end_index; i++) if(col_index[i] <= col) index_col2++;
	

	if(index_row1 < 0) {index_row1++; index_row2++;}
	if(index_row2 > row_end_index) {index_row1--; index_row2--;}
	if(index_col1 < 0) {index_col1++; index_col2++;}
	if(index_col2 > col_end_index) {index_col1--; index_col2--;}
	// #ifdef DEBUG
    //     cout << "index_row1 : " << index_row1 << endl;
    //     cout << "index_row2 : " << index_row2 << endl;
    //     cout << "index_col1 : " << index_col1 << endl;
    //     cout << "index_col2 : " << index_col2 << endl;
    //     cout << "row : " << row << endl;
    //     cout << "col : " << col << endl;
    // #endif
    // xxx_yyy: distance: yyy - xxx
	double col2_col = col_index[index_col2] - col;
	double col2_col1 = col_index[index_col2] - col_index[index_col1];
	double col_col1 = col - col_index[index_col1];
	double row2_row = row_index[index_row2] - row;
	double row2_row1 = row_index[index_row2] - row_index[index_row1];
	double row_row1 = row - row_index[index_row1];
    
	return col2_col * row2_row / col2_col1 / row2_row1 * table[index_row1][index_col1]
		+ col_col1 * row2_row / col2_col1 / row2_row1 * table[index_row1][index_col2]
		+ col2_col * row_row1 / col2_col1 / row2_row1 * table[index_row2][index_col1]
		+ col_col1 * row_row1 / col2_col1 / row2_row1 * table[index_row2][index_col2];
}

void Module::max_delay_path_calculation(int cur_gate_index, double total_output_net_cap, CellLibrary &cell){
	// int max_path_prev_gate_index = -1;
	// double prev_max_delay = 0;
	// double prev_max_transition = 0;
    vector <double> prev_max_delay,prev_max_transition;
    bool Port1,Port2,Input_Port;
	for(int i =0; i < gates[cur_gate_index].output.size(); i++){
        prev_max_delay.push_back(0);
        prev_max_transition.push_back(0);
        if(gates[cur_gate_index].prev.size() == gates[cur_gate_index].wires.size()-1){
            prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.front()][i];
            prev_max_transition[i] = gates[gates[cur_gate_index].prev.front()].critical_transition[i];
            max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.front();

            for(list<int>::iterator lp = gates[cur_gate_index].prev.begin(); lp != gates[cur_gate_index].prev.end(); lp++){
                if(gates[*lp].critical_transition[i] > prev_max_transition[i]){
                    prev_max_transition[i] = gates[*lp].critical_transition[i];
                }
                if(*lp == gates[cur_gate_index].prev.front()){
                    Port1 = gates[*lp].output[i];
                }
                else{
                    Port2 = gates[*lp].output[i];
                }
                // if(max_delay[*lp][i] > prev_max_delay[i]){
                //     prev_max_delay[i] = max_delay[*lp][i];
                //     max_delay_path_from_prev_gate[cur_gate_index][i] = *lp;
                // }
            }

            if(gates[cur_gate_index].type == "NOR2X1"){
                if(Port1 == true && Port2 == false){
                    prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.front()][i];
                    max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.front();
                }
                else if (Port1 == false && Port2 == true ){
                    prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.back()][i];
                    max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.back();
                }
                else if (Port1 == true && Port2 == true ){
                    if(prev_max_delay[i] < max_delay[gates[cur_gate_index].prev.back()][i]){
                        prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.front()][i];
                        max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.front();
                    }
                    else{
                        prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.back()][i];
                        max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.back();
                    }
                    
                }
                else{
                    if(max_delay[gates[cur_gate_index].prev.front()][i] > prev_max_delay[i]){
                        prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.front()][i];
                        max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.front();
                    }
                    else{
                        prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.back()][i];
                        max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.back();
                    }
                }
            }
            else if(gates[cur_gate_index].type == "NANDX1"){
                cout << "--------------" << endl;
                cout << "Gates : " << gates[cur_gate_index].name << endl;
                cout << "max_delay[gates[cur_gate_index].prev.front()][i] :" << max_delay[gates[cur_gate_index].prev.front()][i] << endl;
                cout << "prev_max_delay[i] :" << prev_max_delay[i] << endl;
                cout << "Port1 :" << Port1 << endl;
                cout << "Port2 :" << Port2 << endl;
                cout << "-----------" << endl;
                if(Port1 == false && Port2 == true){
                    prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.front()][i];
                    max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.front();
                }
                else if (Port1 == true && Port2 == false){
                    prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.back()][i];
                    max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.back();
                }
                else if (Port1 == false && Port2 == false ){
                    if(prev_max_delay[i] < max_delay[gates[cur_gate_index].prev.back()][i]){
                        prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.front()][i];
                        max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.front();
                    }
                    else{
                        prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.back()][i];
                        max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.back();
                    }   
                }
                else{
                    if(max_delay[gates[cur_gate_index].prev.front()][i] > prev_max_delay[i]){
                        prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.front()][i];
                        max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.front();
                    }
                    else{
                        prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.back()][i];
                        max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.back();
                    }
                }
            }
        }
        else if(gates[cur_gate_index].prev.size() == gates[cur_gate_index].wires.size()-2 && gates[cur_gate_index].prev.size()!= 0 ){
            prev_max_transition[i] = gates[gates[cur_gate_index].prev.front()].critical_transition[i];
            max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.front();
            for(list<int>::iterator lp = gates[cur_gate_index].prev.begin(); lp != gates[cur_gate_index].prev.end(); lp++){
                if(gates[*lp].critical_transition[i] > prev_max_transition[i]){
                    prev_max_transition[i] = gates[*lp].critical_transition[i];
                }
            }
            // cout << "fuck " << endl;
            for(int j =0 ; j < gates[cur_gate_index].wires.size()-1; j++){
                for (int x = 0; x < input.size(); x++){
                    if(gates[cur_gate_index].wires[j] == input[x]){
                        // cout << gates[cur_gate_index].wires[j] << endl;
                        // cout << input[x] << endl;
                        Input_Port = gates[cur_gate_index].wires_input[i][j];
                        break;
                    }
                }
                break;
            }
            Port2 = gates[gates[cur_gate_index].prev.front()].output[i];
            
            if(gates[cur_gate_index].type == "NOR2X1"){
                if(Input_Port == true && Port2 == false){
                    prev_max_delay[i] = 0;
                    max_delay_path_from_prev_gate[cur_gate_index][i] = -1;
                }
                else if (Input_Port == false && Port2 == true ){
                    prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.back()][i];
                    max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.back();
                }
                else if (Input_Port == true && Port2 == true ){
                    prev_max_delay[i] = 0;
                    max_delay_path_from_prev_gate[cur_gate_index][i] = -1;
                }
                else{
                    prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.front()][i];
                    max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.front();
                }
            }
            else if(gates[cur_gate_index].type == "NANDX1"){
                if(Input_Port == false && Port2 == true){
                    prev_max_delay[i] = 0;
                    max_delay_path_from_prev_gate[cur_gate_index][i] = -1;
                }
                else if (Input_Port == true && Port2 == false){
                    prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.back()][i];
                    max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.back();
                }
                else if (Input_Port == false && Port2 == false ){
                    prev_max_delay[i] = 0;
                    max_delay_path_from_prev_gate[cur_gate_index][i] = -1; 
                }
                else{
                    prev_max_delay[i] = max_delay[gates[cur_gate_index].prev.front()][i];
                    max_delay_path_from_prev_gate[cur_gate_index][i] = gates[cur_gate_index].prev.front();
                }
            }
        }
        // max_path_prev_gate_index = max_delay_path_from_prev_gate[cur_gate_index];
        double max_input_transition = prev_max_transition[i];
        // cout << "max input transistion : " << max_input_transition << endl;
        if (gates[cur_gate_index].output[i]==true){
            double rise_transition = lookup_table(max_input_transition, total_output_net_cap, cell.index_2, cell.index_1, cell.rise_transition);
            double cell_rise = lookup_table(max_input_transition, total_output_net_cap, cell.index_2, cell.index_1, cell.cell_rise);
            gates[cur_gate_index].critical_transition.push_back(rise_transition);
            gates[cur_gate_index].critical_cell_delay.push_back(cell_rise);
            max_delay[cur_gate_index][i] = prev_max_delay[i] + gates[cur_gate_index].critical_cell_delay[i];
	        // cout << "--------------------------" << endl;
            // cout << "Gates : " << gates[cur_gate_index].name << endl;
            // cout << "rise_transition: " << rise_transition << endl;
            // cout << "cell_rise: " << cell_rise << endl;
            // cout << "--------------------------" << endl;
        }
        else{
            double cell_fall = lookup_table(max_input_transition, total_output_net_cap, cell.index_2, cell.index_1, cell.cell_fall);
            double fall_transition = lookup_table(max_input_transition, total_output_net_cap, cell.index_2, cell.index_1, cell.fall_transition);
            gates[cur_gate_index].critical_transition.push_back(fall_transition);
            gates[cur_gate_index].critical_cell_delay.push_back(cell_fall);
            max_delay[cur_gate_index][i] = prev_max_delay[i] + gates[cur_gate_index].critical_cell_delay[i];
            // cout << "--------------------------" << endl;
            // cout << "Gates : " << gates[cur_gate_index].name << endl;
            // cout << "fall_transition: " << fall_transition << endl;
            // cout << "cell_fall: " << cell_fall << endl;
            // cout << "--------------------------" << endl;
        }
    }
	
}

void Module::delay(map<string, CellLibrary> &cells){
	vector <double> tmp_v;
    vector <int> tmp_v2;
    for(int i = 0; i < gates[1].output.size(); i++){
        tmp_v.push_back(0);
        tmp_v2.push_back(-1);
    }
    for(int i = 0; i < gates.size(); i++){
		max_delay.push_back(tmp_v);
		max_delay_path_from_prev_gate.push_back(tmp_v2);
	}
    // caculate output capatiance
	for(int i = 0; i < topological_order.size(); i++){
		int cur_gate_index = topological_order[i];
		double total_output_net_cap = (gates[cur_gate_index].isOutput) ? 0.03 : 0.0;
		for(list<int>::iterator lp = gates[cur_gate_index].next.begin(); lp != gates[cur_gate_index].next.end(); lp++){
			for(int k = 0; k < gates[*lp].wires.size()-1; k++){
                //current output cap == next input cap
                if(gates[cur_gate_index].wires.back() == gates[*lp].wires[k]){
					total_output_net_cap += cells[gates[*lp].type].pins[gates[*lp].ports[k]].cap;
					break;
				}
			}
		}
        // #ifdef DEBUG
        //     cout << "--------------------------" << endl;
        //     cout << "topological_order : "<< cur_gate_index << endl;
        //     cout << "gate : " << gates[cur_gate_index].name << endl;
        //     cout << "Out cap :" << total_output_net_cap << endl;
        // #endif
		max_delay_path_calculation(cur_gate_index, total_output_net_cap, cells[gates[cur_gate_index].type]);
        // cout << "--------------------------" << endl;
	}
	
	critical_path_out_index = output_gates.front();
    for(int x = 0; x < input_signal.size(); x++){ // patten number
        for(int i = 1; i < output_gates.size(); i++){
            if(max_delay[output_gates[i]][x] > max_delay[critical_path_out_index][x])
                critical_path_out_index = output_gates[i];
        }
    }
}

void Module::output_file(){
	fstream fout;
	fout.open("Out_ans.txt", ios::out);
    
    for(int j=0; j < input_signal.size(); j++){
        fout << "Longest delay = " << max_delay[critical_path_out_index][j] << ", the path is:" << endl;
        vector<string> max_delay_path;
        int index = critical_path_out_index;
        while(max_delay_path_from_prev_gate[index][j] != -1){
        	max_delay_path.push_back(gates[index].wires.back());
        	index = max_delay_path_from_prev_gate[index][j];
        }
        max_delay_path.push_back(gates[index].wires.back());
        max_delay_path.push_back(gates[index].wires.front());
        for(int i = max_delay_path.size()-1; i > 0; --i)
        	fout << max_delay_path[i] << " -> ";
        fout << max_delay_path.front() << endl;
        fout << endl;
        for (int i=0; i < gates.size(); i++){
            fout << gates[i].name << " "
                << gates[i].output[j] << " "
                << gates[i].critical_cell_delay[j] << " "
                << gates[i].critical_transition[j] << endl;
        }
        fout << endl;
    }
	fout.close();
}


class Parser {
public:
    void parse_lu_table(fstream &infile, vector<double> &idx1, vector<double> &idx2) {
        string value, tmp;
        int begin, end;
        stringstream ss;
        while(getline(infile, tmp, ';')){
            if(tmp.find("index_1") != string::npos){
                begin = tmp.find('"');
                end = tmp.find_last_of('"');
                value = tmp.substr(begin + 1, end - begin - 1);
                ss.str("");
                ss.clear();
                ss << value;
                while (getline(ss, tmp, ',')) {
                    idx1.push_back(atof(tmp.c_str()));
                }
                ss.str("");
                ss.clear();
            }
            if(tmp.find("index_2") != string::npos){
                begin = tmp.find('"');
                end = tmp.find_last_of('"');
                value = tmp.substr(begin + 1, end - begin - 1);
                ss.str("");
                ss.clear();
                ss << value;
                while (getline(ss, tmp, ',')) {
                    idx2.push_back(atof(tmp.c_str()));
                }
                ss.str("");
                ss.clear();
                break;
            }
        }
    }

    void parse_cell(fstream &fin, CellLibrary &c) {
        string in, temp;
        int num;
        stringstream ss1, ss2;
        int cnt = 0;
        vector<double> value;
        while (cnt < 4) {
            getline(fin, in);
            ss2 << in;
            // cout << "in/ss2 :" << in << endl;
            ss2 >> temp;
            // cout << "before tmp/ss1 :" << temp << endl;
            ss1 << temp;
            getline(ss1, temp, '(');
            // cout << "After tmp :" << temp << endl;
            if (temp == "pin") {
                Pin p;
                getline(ss1, p.name, ')');
                // cout << "pname :" << ss2.str() << endl;
                parse_pin(fin, p);
                // cout << "pname :" << p << endl;
                c.pins[p.name] = p;
            } else if (temp == "cell_rise") {
                parse_timing(fin, c, temp);
                cnt++;
            } else if (temp == "cell_fall") {
                parse_timing(fin, c, temp);
                cnt++;
            } else if (temp == "rise_transition") {
                parse_timing(fin, c, temp);
                cnt++;
            } else if (temp == "fall_transition") {
                parse_timing(fin, c, temp);
                cnt++;
            }

            ss1.str("");
            ss1.clear();
            ss2.str("");
            ss2.clear();
        }
        // #ifdef DEBUG
        //     cout << "--------------------------" << endl;
        //     cout << "cell rise" << endl;
        //     cout << endl;
        //     vector< vector<double> >::iterator item2v = c.cell_rise.begin();
        //     while(item2v !=c.cell_rise.end()){
        //         vector<double>::iterator item = (*item2v).begin();
        //         while(item != (*item2v).end()){
        //             cout<<*item<<"\t";
        //             item++;
        //         }
        //         cout << endl;
        //         item2v++;
        //     }
        //     cout << "--------------------------" << endl;
        // #endif
            
    }

    void parse_pin(fstream &fin, Pin &p) {
        string in, temp;
        int begin, end;
        getline(fin, temp);
        begin = temp.find(':');
        end = temp.find_last_of(';');
        in = temp.substr(begin + 2, end - (begin + 2));
        p.direct = in;
        getline(fin, temp);
        begin = temp.find(':');
        end = temp.find_last_of(';');
        in = temp.substr(begin + 2, end - (begin + 2));
        p.cap = atof(in.c_str());
    }
        
    void parse_timing(fstream &fin, CellLibrary &c, string timing) {
        string in, tmp;
        int begin, end;
        stringstream ss;
        
        for (int i = 0; i < c.index_1.size(); i++) {
            getline(fin, tmp);
            begin = tmp.find('"');
            end = tmp.find_last_of('"');
            in = tmp.substr(begin + 1, end - begin - 1);
            ss << in;
            vector<double> value;
            while (getline(ss, tmp, ',')) {
                value.push_back(atof(tmp.c_str()));
            }
            if (timing == "cell_rise"){
                c.cell_rise.push_back(value);
            }
            else if (timing == "cell_fall"){
                c.cell_fall.push_back(value);
            }
            else if (timing == "rise_transition"){
                c.rise_transition.push_back(value);
            }
            else if (timing == "fall_transition"){
                c.fall_transition.push_back(value);
            }
            ss.str("");
            ss.clear();
        }
    }
    
    void parse_Lib(char *file, map<string, CellLibrary> &cells) {
        fstream infile;
        infile.open(file, ios::in);
        if (!infile) {
            cerr << "Fail to open: " << file << endl;
            exit(1);
        }

        string in, tmp;
        int num;
        stringstream ss;
        vector<double> output_cap, input_tran;

        while (getline(infile, in)) {
            ss << in;
            while (getline(ss, tmp, '(')) {
                if (tmp == "lu_table_template") {
                    parse_lu_table(infile, output_cap, input_tran);
                } else if (tmp == "cell ") {
                    CellLibrary c;
                    c.index_1 = output_cap;
                    c.index_2 = input_tran;
                    getline(ss, c.name, ')');
                    parse_cell(infile, c);
                    cells[c.name] = c;
                }
            }

            ss.str("");
            ss.clear();
        }

        infile.close();
    }
    
    void parse_pat(char *file, Module &m) {
        fstream infile;
        infile.open(file, ios::in);
        if (!infile) {
            cerr << "Fail to open: " << file << endl;
            exit(1);
        }
        string tmp;
        stringstream ss;
        vector<string> input_name;
        vector<bool> input_signal;
        vector<bool> tmp_v;
        
        // input name
        getline(infile, tmp);
        ss << tmp;
        ss >> tmp;
        if(tmp == "input"){
            while (getline(ss, tmp, ',')){
                tmp.erase(0,1);
                input_name.push_back(tmp);
            }
        }
        ss.str("");
        ss.clear();	

        // input signal
        while(getline(infile, tmp)){
            ss << tmp;
            while (ss >> tmp){
                if (tmp == ".end"){
                    if(input_signal.size() == input_name.size()){
                        for(int i=0 ; i < m.input.size() ;i++){
                            for(int j=0 ; j < input_name.size() ;j++){
                                if(input_name[j] == m.input[i]){
                                    tmp_v.push_back(input_signal[j]);
                                    break;
                                }
                            }
                        }
                        m.input_signal.push_back(tmp_v);
                        input_signal.clear();
                        tmp_v.clear();
                    }
                    break;
                }
                else{
                    if(input_signal.size() == input_name.size()){
                        for(int i=0 ; i < m.input.size() ;i++){
                            for(int j=0 ; j < input_name.size() ;j++){
                                if(input_name[j] == m.input[i]){
                                    tmp_v.push_back(input_signal[j]);
                                    break;
                                }
                            }
                        }
                        m.input_signal.push_back(tmp_v);
                        input_signal.clear();
                        tmp_v.clear();
                        if(tmp == "1"){
                            input_signal.push_back(true);
                        }
                        else{
                            input_signal.push_back(false);
                        }
                    }
                    else{
                        if(tmp == "1"){
                            input_signal.push_back(true);
                        }
                        else{
                            input_signal.push_back(false);
                        }
                    }
                }
            }
            ss.str("");
            ss.clear();	
        }
        
        #ifdef DEBUG
            cout << "--------------------------" << endl;
            cout << "input patten order" << endl;
            for(int i=0 ; i < m.input.size() ;i++){
                cout << m.input[i] << "\t";
            }
            cout << endl;
            vector< vector<bool> >::iterator item2v = m.input_signal.begin();
            while(item2v !=m.input_signal.end()){
                vector<bool>::iterator item = (*item2v).begin();
                while(item != (*item2v).end()){
                    cout<<*item<<"\t";
                    item++;
                }
                cout << endl;
                item2v++;
            }
            cout << "--------------------------" << endl;
        #endif
        
    }

    void parse_verilog(char *file, Module &m) {
        ifstream infile;
        infile.open(file, ios::in);
        if (!infile) {
            cerr << "Fail to open: " << file << endl;
            exit(1);
        }

        string tmp;
        stringstream ss;
        // module
        while(getline(infile, tmp, ';')){
            ss << tmp;
            ss >> tmp;
            if(tmp == "module"){
                m.name = tmp;
                ss.str("");
                ss.clear();			
			}
            // cout << "tmp :" << tmp << endl;
            else if (tmp == "input"){
                while (getline(ss, tmp, ',')){
                    tmp.erase(0,tmp.find_first_not_of(" "));
                    m.input.push_back(tmp);
                }
                ss.str("");
	            ss.clear();
            }
            else if (tmp == "output"){
                while (getline(ss, tmp, ',')){
                    tmp.erase(0,tmp.find_first_not_of(" "));
                    m.output.push_back(tmp);
                }
                ss.str("");
	            ss.clear();
            }
            else if (tmp == "wire"){
                while (getline(ss, tmp, ',')){
                    tmp.erase(0,tmp.find_first_not_of(" "));
                    m.wire.push_back(tmp);
                }
                ss.str("");
	            ss.clear();
            }
            else if (tmp == "endmodule"){
                break;
            }
            else{
                Gate g;
                g.type = tmp;
                while (ss >> tmp) {
                    g.name = tmp;
                    while (getline(ss, tmp, ',')) {
                        // cout << "tmp : " << tmp << endl;
                        tmp.erase(0,tmp.find_first_not_of(" "));
                        if (tmp[0] == '.') {
                            int start = tmp.find("(");
                            int end = tmp.find(")");
                            // cout << "tmp :" << tmp.substr(start + 1, end - start - 1) << endl;
                            g.ports.push_back(tmp.substr(1, start - 1));
                            g.wires.push_back(tmp.substr(start + 1, end - start - 1));
                        }
                        else if (tmp[0] == '(') {
                            tmp.erase(0,1);
                            int start = tmp.find("(");
                            int end = tmp.find(")");
                            // cout << "tmp :" << tmp.substr(1, start - 1) << endl;
                            g.ports.push_back(tmp.substr(1, start - 1));
                            g.wires.push_back(tmp.substr(start + 1, end - start - 1));
                        }
                    }
                    m.gates.push_back(g);
                }
                ss.str("");
                ss.clear();
            }
        }
        infile.close();
    }

};

int main(int argc, char *argv[]){
    map<string, CellLibrary> cells;
    Parser parser;
    Module module;
    string test;
    parser.parse_Lib(argv[5], cells);
    parser.parse_verilog(argv[1], module);
    parser.parse_pat(argv[3], module);
    module.topological_sort();
    module.delay(cells);
    module.output_file();
    return 0;
};