# Snakmake

Оу, всем привет. Сегодня мы расмотрим технологию Snakemake (SNK).

Подозреваю у каждого были моменты, когда вы последовательно обрабатывали данные, и вдруг поняли, что ошиблись с командой в самом начале анализа. И так как, каждый следующий шаг зависил от предыдущего, вам приходилось выполнять всю последовательность команд заново. Согласитесь неудобно? Подобную проблему решает SNK,  позволяет создавать пайплайны из разных сниппетов кода, причем можно объединять любые скрипты (Bash, R, Python).

Да, если вы любите красивые IDE, то у Pycharm есть специальный плагин для работы со SNK - [https://plugins.jetbrains.com/plugin/11947-snakecharm](https://plugins.jetbrains.com/plugin/11947-snakecharm) 

SNK преследует парадигму последовательного выполнения команд, пайплайн складывается в виде набора правил (rule) - как определить output файл из input файла. Зависимости между такими правилами определяются автоматически - создается направленный ациклический граф вычислений, которые могут быть параллелизованы.

Давайте установим SNK и попробуем сделать первую команду.

!pip install snakemake 

Сами команды SNK пишутся в отдельном файле Snakefile с указанием числа используемых ядер.

Напишем первую команду записи строки в файл.

```bash
rule generate:
    output: 'words.txt'
    shell: 'echo "Hello, world!" > {output}'
```

```bash
snakemake -s Snakefile --cores all
```

В этом примере создается текстовый файл с именем words.txt .

Давайте добавим правило для подсчета символов в строке.

```bash
rule count_char:
    input: 'words.txt'
    output: 'len_string.txt'
    shell: 'wc -c {input} > {output}'

rule generate:
    output: 'words.txt'
    shell: 'echo "Hello, world!" > {output}'
```

SNK выполняет вычисления экономно, то есть если в директории уже есть один из промежуточных файлов, которые вы получаете во время выполнения пайплайна SNK процесс вычислений не запустится. Для игнорирования этого, стоит добавить флаг -f или удалить все файлы полученные до этого.

Покольку файл 'words.txt' уже есть:

```bash
snakemake -s Snakefile -f --cores all
```

Как задается порядок команд? 

Порядок задаете вы - логика следующая: если в правиле А - output файл это x, а

а правиле B - input файл тоже x, то правило B идет за правилом A.

Пример.

```bash
rule convert_to_bam:
    input:
        "input.sam"
    output:
        "output.bam"
    shell:
        "samtools view -S -b {input} > {output}"

rule generate_pileup:
    input:
        bam="output.bam",
        reference="reference.fasta"
    output:
        "output.pileup"
    shell:
        "samtools mpileup -f {input.reference} {input.bam} > {output}"
```

Сначала получаем bam файл, затем отправляем его для получения 

**Использование Python скриптов в Snakemake.**

Представим задачу - допустим есть файл со строкой, нужно посчитать кол-во уникальных элементов в строке и записать в отдельный файл. Для этого помогут пара функций на Python.

```bash
def calculate_ch(string):
    stat = {}
    for ch in string:
        if ch not in stat:
            stat[ch] = 1
        else:
            stat[ch] += 1

    ans = []
    for ch in stat:
        entry = ch + '-' + str(stat[ch])
        ans.append(entry)
    return ' '.join(ans)

def process_data(input_file, output_file):
    with open(input_file) as inp:
        string = inp.read()
        ans = calculate_ch(string)

    with open(output_file, 'w') as out:
        out.write(ans)

rule process_rule:
    input:
        "input_files/{file}"
    output:
        "output_files/{file}"
    run:
        process_data(input[0], output[0])
```

```bash
snakemake -s snakefile.smk --cores 1
```

Python скрипты запускаются с помощью метода run.

Ок, а что если таких файлов не один, а много?. В этом помогут wildcards.

**Wildcards** полезны, когда у вас есть набор связанных файлов, которые следуют шаблону имени и вы хотите выполнить одно и то же правило для всех них.  Вы можете написать одно правило, которое применяется к нескольким файлам, без явного перечисления каждого файла.

Пусть в директории `input_files` лежат файлы, и мы хотим над каждым сделать операцию подсчета уникальных символов и записать результат в директорию `output_files.`

Все тоже самое, как в предыдщем примере, теперь нужна заготовка с перечислением всех файлов идущих правилу на вход.

Правило all как бы собирает результат всех вычислений, которые мы получили в правиле process_rule. 

```bash
import os
files = os.listdir('input_files')

rule all:
    input:
        [f'output_files/{file}' for file in files]

rule process_rule:
    input:
        "input_files/{file}"
    output:
        "output_files/{file}"
    run:
        process_data(input[0], output[0])
```

Иной способ, как можно склеить все input файлы - это метод expand.

```bash

import os
files = os.listdir('input_files')

rule all:
    input:
        expand("output_files/{file}", file = files)

rule process_rule:
    input:
        "input_files/{file}"
    output:
        "output_files/{file}"
    run:
        process_data(input[0], output[0])
```

**Распараллеливание.**

Да, поскольку SNK параллелит исполнение команд, можно  задавать число ядер как для исполнения всего пайплайна.

Что интересно, у меня не исполняется snakefile если не указывать число ядер на исполение. 

```bash
snakemake -s snakefile.smk
```

*Error: you need to specify the maximum number of CPU cores to be used at the same time. If you want to use N cores, say --cores N or -cN. For all cores on your system (be sure that this is appropriate) use --cores all.* 

Также возможно  задавать число потоков на одно правило.

```bash
rule NAME:
    input: "path/to/inputfile"
    output: "path/to/outputfile"
    threads: 8
    shell: "somecommand --threads {threads} {input} {output}"
```

Давайте решим задачу выравнивания коротких прочтений  на референс с помощью  snakemake.
Input: reads.fastq.gz ( данные прочтений)
Output: reads.mpileup (файл variant calling на на референсный геном)
Snakefile: snakefile.smk

Выделим набор правил:

1. downoload_ref - загрузка референса.
2. index_ref - индексация референса 
3. align - выравнивание коротких прочтений 
4. index_bam - идексация прочтений 
5. get_var - получение вариантов на позицию. 

**Фишка** 

поскольку SNK составляет из набора правил граф отношений - его можно отрисовать. 

Для это нужна библиотека - *graphviz.*

```bash
snakemake  -s snakefile.smk --cores all -f --dag | dot -Tsvg > pipeline.svg
```

![My Image](https://github.com/GavrilenkoA/snakemake/blob/main/pipeline.svg)
В итоге весь набор правил отрисуется в виде графа - пригодится, если вам в презентации нужно будет показать пайплайн вашего анализа, будет альтернативой выковыриванию стролочек в Power Point.

Я рассмотрел только небольшую функций SNK, для большего ознакомления советую почитать документацию, она довольно приятно написана.

Источники:

1. Курс ****Управление вычислениями -**** [https://stepik.org/course/1612/info](https://stepik.org/course/1612/info)
2. **Документация Snakemake** - [https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html)
