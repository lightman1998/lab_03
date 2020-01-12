#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

class SimplexMethod {
	// строка таблицы
	struct Row {
		vector<double> a; // коэффициенты при x
		double b; // правая часть
	};

	int n; // число переменных
	int m; // число ограничений
	vector<Row> table; // симплекс таблица
	vector<double> c; // коэффициенты оптимизируемой функции
	vector<int> variables; // все переменные
	vector<int> basis; // базисные переменные
	vector<double> deltas; // дельты
	vector<double> solvation;

	void CalculateDeltas(); // вычисление дельт

	int GetArgMinDelta(); // вычисление номера минимальной дельты
	int GetArgMaxDelta(); // вычисление номера максимальной дельты

	void InitialVariables(); // инициализация переменных и базиса

	void MakeNewBasis(double pivot, int index, int jindex); // задание нового базисного элемента
	int MaxNegativeB(); // максимальная по модулю отрицательная b

	int CheckForIntegerSolvation(); // проверка на то, что решение целое
	void CheckNewCondition(vector<int> &best, int &f, int basisIndex); // проверка введённых условий
	bool CheckRestriction(vector<int> solve); // проверка ограничений для решения
	void GetIntegerSolves(vector<vector<int>> &solves, vector<int> &f, vector<int> &solve, int index); // получение полного набора решений(рекурсивно)
	int GetF(vector<int> solve); // получение значения целевой функции по вектору
	void PrintVector(vector<int> vec); // вывод вектора
	void Check(vector<int> &best, int &bestF, int basisIndex);

	bool isSolve;
public:
	SimplexMethod(int index = -1, int resist = 0); // конструктор по умолчанию со всеми данными
	SimplexMethod(const SimplexMethod& method); // конструктор для двойственной задачи из прямой

	void Read(); // ввод значений
	void Print(); // вывод таблицы

	void Solve(int max); // решение ОЗЛП
	void RemoveNegativeB(); // удаление отрицательных b

	void IntegerSolve(); // поиск целочисленного решения
};

// конструктор по умолчанию со всеми данными
SimplexMethod::SimplexMethod(int index, int resist) {
	// задаём значения для конкрутной задачи
	n = 3;
	m = 3 + (index >= 0 ? 1: 0); // увеличиваем кол-во условий, если нужно

	table = vector<Row>(m, { vector<double>(n), 0 });
	c = vector<double>(n);

	c = {2, 6, 7};

	table[0].a = {3, 1, 1};
	table[0].b = {3};

	table[1].a = {1, 2, 0};
	table[1].b = {8};
	
	table[2].a = {0, 0.5, 2};
	table[2].b = {1};
	solvation = vector<double>(n);

	// заполняем новое условие
	if (index >= 0) {
		for (int i = 0; i < n; i++) 
			if (index == i)
				table[m - 1].a[i] = resist < 0 ? -1 : 1;
			else
				table[m - 1].a[i] = 0;

		table[m - 1].b = resist;
	}

	InitialVariables(); // инициализируем переменные
	RemoveNegativeB();

	isSolve = false;
}

// конструктор для двойственной задачи из прямой
SimplexMethod::SimplexMethod(const SimplexMethod& method) {
	// запоминаем размеры
	n = method.n;
	m = method.m;

	table = vector<Row>(m, { vector<double>(n), 0 }); // создаём таблицу
	c = vector<double>(n); // создаём вектор коэффициентов

	// функции цели присваиваем значения правой части прямой задачи
	for (int i = 0; i < n; i++) 
		c[i] = method.table[i].b;

	// правой части присваиваем значения функции цели прямой задачи * -1(так как сразу меняем знак для приведения к неравенствам вида <=)
	for (int i = 0; i < m; i++)
		table[i].b = method.c[i] * -1;

	// транспонируем матрицу прямой задачи и умножаем элементы на -1 (так как сразу меняем знак для приведения к неравенствам вида <=)
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			table[i].a[j] = method.table[j].a[i] * -1;
		}
	}

	solvation = vector<double>(3);
	InitialVariables(); // инициализируем переменные
}


// инициализация переменных и базиса
void SimplexMethod::InitialVariables() {
	variables.clear(); // очищаем переменные

	// добавляем переменные
	for (int i = 0; i < n; i++)
		variables.push_back(i);

	for (int i = 0; i < m; i++) {
		c.push_back(0); // добавляем нули в функцию
		variables.push_back(n + i); // добавляем доп переменные
		basis.push_back(n + i); // делаем их базисными

		// добавляем коэффициенты для переменных с коэффициентом 1, если они стоят на главной диагонали, иначе с нулём
		for (int j = 0; j < m; j++)
			table[i].a.push_back(i == j);
	}
}


// поиск макисмальной отрицательной b
int SimplexMethod::MaxNegativeB() {
	int imax = -1;

	// ищем максимальный отрицательный элемент
	for (int i = 1; i < m; i++) {
		if (table[i].b < 0 && (imax == -1 || table[i].b < table[imax].b))
			imax = i;
	}

	return imax; // возвращаем максимум
}

// устранение отрицательной правой части
void SimplexMethod::RemoveNegativeB() {

	int imax = MaxNegativeB(); // индекс максимального по модулю отрицательного элемента

	// пока если отрицательные элементы
	while (imax != -1) {
		int jmax = 0; // индекс новой базисной переменной

		// идём по столбцу и ищем максимальный по модул. элемент
		for (int j = 1; j < m; j++) {
			if (fabs(table[imax].a[j]) > fabs(table[imax].a[jmax]))
				jmax = j;
		}

		basis[imax] = jmax;	// запоминаем индекс новой базисной переменной
		MakeNewBasis(table[imax].a[jmax], imax, jmax); // делаем этот элемент базисным

		imax = MaxNegativeB(); // находим новый максимальный по модулю элемент в правой части
	}
}

// создание новой базисной переменной на месте index, jindex
void SimplexMethod::MakeNewBasis(double pivot, int index, int jindex) {
	// делим строку на элемент
	for (size_t i = 0; i < table[index].a.size(); i++)
		table[index].a[i] /= pivot;

	table[index].b /= pivot;
	
	// вычитаем из всех остальных строк эту строку, умноженную на элемент в столбце jmax
	for (int i = 0; i < m; i++) {
		if (i == index)
			continue;

		double value = table[i].a[jindex];

		for (size_t j = 0; j < table[i].a.size(); j++)
			table[i].a[j] -= table[index].a[j] * value;

		table[i].b -= table[index].b * value;
	}
}


// ввод значений
void SimplexMethod::Read() {
	cout << "Enter function coefficients (c): ";
	c = vector<double>(n); // создаём вектор коэффициентов

	// считываем коэффициенты оптимизируемой функции
	for (int i = 0; i < n; i++)
		cin >> c[i];

	cout << "Enter restrictions coefficients:" << endl;

	// считываем коэффициенты ограничений
	for (int i = 0; i < m; i++) {
		cout << "Enter restriction " << (i + 1) << ": ";

		for (int j = 0; j < n; j++)
			cin >> table[i].a[j];

		cin >> table[i].b;
	}

	variables.clear(); // очищаем переменные

	// добавляем переменные
	for (int i = 0; i < n; i++)
		variables.push_back(i);

	for (int i = 0; i < m; i++) {
		c.push_back(0); // добавляем нули в функцию
		variables.push_back(n + i); // добавляем доп переменные
		basis.push_back(n + i); // делаем их базисными

		// добавляем коэффициенты для переменных с коэффициентом 1, если они стоят на главной диагонали, иначе с нулём
		for (int j = 0; j < m; j++)
			table[i].a.push_back(i == j);
	}
}

// вывод таблицы
void SimplexMethod::Print() {
	int vars = variables.size();

	cout << endl;
	cout << "+-----+";

	for (int i = 0; i < vars; i++)
		cout << "-----------+";

	cout << endl;

	cout << "|  C  |";

	for (int i = 0; i < vars; i++)
		cout << " " << setw(9) << c[i] << " |";

	cout << endl;

	cout << "+-----+";

	for (int i = 0; i <= vars; i++)
		cout << "-----------+";

	cout << endl;

	cout << "|basis|";
	for (int i = 0; i < vars; i++)
		cout << "    x" << setw(2) << left << (i + 1) << "    |";

	cout << "     b     |" << endl;
	cout << "+-----+";

	for (int i = 0; i <= vars; i++)
		cout << "-----------+";

	cout << endl;

	for (int i = 0; i < m; i++) {
		cout << "| x" << setw(2) << left;

		if ((size_t)i < basis.size())
			cout << (basis[i] + 1);
		else
			cout << "?";

		cout  << " |";

		for (size_t j = 0; j < table[i].a.size(); j++)
			cout << " " << setw(9) << table[i].a[j] << " |";

		cout << " " << setw(9) << table[i].b << " |" << endl;
	}

	cout << "+-----+";

	for (int i = 0; i <= vars; i++)
		cout << "-----------+";

	cout << endl;

	if (!deltas.size())
		return;

	cout << "|  D  |";

	for (size_t i = 0; i < deltas.size(); i++)
		cout << " " << setw(9) << deltas[i] << " |";

	cout << endl;

	cout << "+-----+";

	for (int i = 0; i <= vars; i++)
		cout << "-----------+";

	cout << endl;


	cout << endl << endl << "X = (";

	for (int i = 0; i < n; i++) {
		bool inBasis = false;

		for (int j = 0; j < m; j++) {
			if (basis[j] == i) {
				cout << table[j].b << " ";
				inBasis = true;
				solvation[i] = table[j].b;
				break;
			}
		}

		if (!inBasis) {
			solvation[i] = 0;
			cout << "0 ";
		}
	}

	cout << ")";

	cout << "\tF = " << deltas[n + m] << endl << endl;
}

// вычисление дельт
void SimplexMethod::CalculateDeltas() {
	deltas.clear(); // очищаем массив дельт

	// проходимся по всем переменным
	for (size_t i = 0; i <= variables.size(); i++) {
		double delta = 0;

		// вычилсяем дельту
		for (size_t j = 0; j < basis.size(); j++)
			delta += c[basis[j]] * (i < variables.size() ? table[j].a[i] : table[j].b);

		// вычитаем коэффициент функции
		if (i < variables.size())
			delta -= c[i];

		deltas.push_back(delta); // добавляем дельту в массив
	}
}

// вычисление номера минимальной дельты
int SimplexMethod::GetArgMaxDelta() {
	int imax = 0; // считаем, что первая дельта максимальна

	// проходимся по всем дельтам
	for (size_t i = 1; i < deltas.size() - 1; i++)
		if (deltas[i] > deltas[imax]) // если дельта стала больше максимальной
			imax = i; // обновляем индекс максимума

	return imax; // возвращаем индекс максимума
}

// вычисление номера минимальной дельты
int SimplexMethod::GetArgMinDelta() {
	int imin = 0; // считаем, что первая дельта минимальная

	// проходимся по всем дельтам
	for (size_t i = 1; i < deltas.size() - 1; i++)
		if (deltas[i] < deltas[imin]) // если дельта стала меньше минимальной
			imin = i; // обновляем индекс минимума

	return imin; // возвращаем индекс минимума
}

// решение ОЗЛП  max = 1 при минимизации max = -1 при максимизации
void SimplexMethod::Solve(int max) {
	while (true) {
		CalculateDeltas(); // рассчитываем дельты
		int jmax;

		// если минимизируем
		if (max == 1)
			jmax = GetArgMaxDelta(); // ищем индекс максимальной
		// если максимизация
		else
			jmax = GetArgMinDelta(); // ищем индекс минимальной

		double maxDelta = deltas[jmax]; // получаем максимальную дельту

		// если она не положительна(или неотрицатльная для максимизации)
		if (maxDelta * max <= 0) {
			Print(); // выводим таблицу 
			isSolve = true;
			break; // и выходим
		}

		vector<double> Q(m); // создаём симплекс отношения
		int imin = -1;

		// идём по ограничениям
		for (int i = 0; i < m; i++) {
			if (table[i].a[jmax] == 0) { // если коэффициент равен 0
				Q[i] = 0; // то отношение равно нулю
			}
			else {
				Q[i] = table[i].b / table[i].a[jmax]; // вычисляем результат отношения

				// если оно отрицательно, то идём дальше
				if (Q[i] < 0)
					continue;

				// иначе обновляем минимальное симплекс отношение
				if (imin == -1 || Q[i] < Q[imin])
					imin = i;
			}
		}

		if (imin == -1) {
			cout << "Solution does not exist" << endl;
			isSolve = false;
			break;
		}

		basis[imin] = jmax; // делаем переменную базисноц
		double pivot = table[imin].a[jmax]; // получаем опорный элемент
		
		MakeNewBasis(pivot, imin, jmax);
	}	
}

// проверка ограничений для решения
bool SimplexMethod::CheckRestriction(vector<int> solve) {
	// идём по ограничениям
	for (int i = 0; i < m; i++) {
		double sum = 0; // сумма

		// находим сумму
		for (int j = 0; j < n; j++)
			sum += table[i].a[j] * solve[j];

		// если не удовлетворяет ограничению
		if (sum > table[i].b)
			return false; // то не подходит решение
	}

	return true;
}

// получение значения целевой функции по вектору
int SimplexMethod::GetF(vector<int> solve) {
	int res = 0;

	for (int i = 0; i < n; i++)
		res += c[i] * solve[i];

	return res;
}

// получение полного набора решений(рекурсивно)
void SimplexMethod::GetIntegerSolves(vector<vector<int>> &solves, vector<int> &f, vector<int> &solve, int index) {
	// если дошли до конца
	if (index == n) {
		// проверяем ограничения
		if (CheckRestriction(solve)) {
			solves.push_back(solve); // закидываем в решения
			f.push_back(GetF(solve)); // получаем значение функции цели
		}

		return; // выходим
	}

	// рекурсивно запускаем программу для каждой следующей переменной
	for (int i = 0; i < 100; i++) {
		solve[index] = i;

		GetIntegerSolves(solves, f, solve, index + 1);
	}
}

// вывод вектора
void SimplexMethod::PrintVector(vector<int> vec) {
	cout << "(";
		
	for (size_t i = 0; i < vec.size() - 1; i++) 
		cout << vec[i] << ",";

	cout << vec[vec.size() - 1] << ")";
}

// проверка на то, что решение целое
int SimplexMethod::CheckForIntegerSolvation() {
	if (!isSolve)
		return -1;

	for (int i = 0; i < n; i++) 
		if (solvation[i] - floor(solvation[i]) != 0)
			return 0;

	return 1;
}

// проверка введённых условий
void SimplexMethod::CheckNewCondition(vector<int> &best, int &bestF, int basisIndex) {
	cout << "Initial new simplex table" << endl;
	Print(); // вывод исходной симплекс-таблицы

	cout << "Result table:" << endl;
	Solve(-1); // решение(там и вывод результата)

	int checkResult = CheckForIntegerSolvation();
	// если решение целое
	if (checkResult == 1) {
		cout << "Since the answer is only integers, we consider this a solution" << endl; // выводим, что всё норм
	
		// если целевая функция больше лучшей
		if (deltas[n + m] > bestF) {
			bestF = deltas[n + m]; // запоминаем новое лучшее значение
			
			// заполняем вектор
			for (int i = 0; i < n; i++)
				best[i] = (int)solvation[i];
		}
	}
	// иначе запускаем для других переменных
	else if (checkResult == 0) {
		cout << "The answer is not only integers" << endl;
		Check(best, bestF, basisIndex);
	}
}

// рекурсивная проверка всех ветвей дерева
void SimplexMethod::Check(vector<int> &best, int &bestF, int basisIndex) {
	// идём по всем переменным
	for (int i = 0; i < n; i++) {
		if (i == basisIndex)
			continue;

		// если переменная целая, то пропускаем её
		if (solvation[i] - floor(solvation[i]) == 0) {
			cout << "Variable x" << (i + 1) << " is integer => it is solvation" << endl;
			continue;
		}

		cout << "Variable x" << (i + 1) << " isn`t integer. Consider 2 cases: " << endl;
		cout << "First: x" << (i + 1) << "<= " << floor(solvation[i]) << endl << endl; // вывод первого ветвления

		SimplexMethod method1(i, floor(solvation[i])); // создание таблицы с новым ограничением
		method1.CheckNewCondition(best, bestF, i); // решение для таких ограничений

		cout << endl << endl << "Second: x" << (i + 1) << ">= " << (floor(solvation[i]) + 1) << endl << endl; // вывод 2 ветвления

		SimplexMethod method2(i, -(floor(solvation[i]) + 1)); // создание таблицы с новым ограничением
		method2.CheckNewCondition(best, bestF, i); // решение для таких ограничений

		cout << endl << endl;
	}
}

// поиск целочисленного решения
void SimplexMethod::IntegerSolve() {
	vector<vector<int>> solves; // полнй набор решений

	vector<int> solve(n, 0); // решение
	vector<int> f; // вектор функций цели для всех решений

	GetIntegerSolves(solves, f, solve, 0); // получение всех решений

	cout << "All integer solves:" << endl;

	int imax = 0; // индекс лучшего

	// идём по наборам
	for (size_t i = 0; i < solves.size(); i++) {
		PrintVector(solves[i]); // выводим

		// если нашли лучше максимума
		if (f[i] > f[imax])
			imax = i; // меняем максимум

		cout << "\tF = " << f[i] << endl; // вывод функции цели
	}

	cout << "Solve with max cost function: ";

	PrintVector(solves[imax]); // вывод лучшего
	cout << "\tF = " << f[imax] << endl << endl; // и его целевой функции

	// поиск решения(не целочисленного)
	cout << "Find solvation by simplex method" << endl;
	cout << "Initial table:" << endl;
	Print(); // вывод исходной таблицы

	cout << "Result table:" << endl;
	Solve(-1); // рашение и вывод результата

	vector<int> best(n); // лучший набор
	int bestF = 0; // лучшая целевая функция

	Check(best, bestF, basis.size()); // проверка всех ветвей решения

	cout << "*********************************RESULT********************************************" << endl;

	cout << endl << "Best integer solvation: ";
	PrintVector(best); // вывод лучшего

	cout << "\tF = " << GetF(best) << endl << endl; // и его функции цели
}

int main() {
	SimplexMethod method; // прямой метод
	method.IntegerSolve();
}
