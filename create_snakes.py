import os
import sys

if __name__ == '__main__':
    inputs = sys.argv[1:]
    with open('Snakefile_template', 'r') as f:
        snake_template = f.read()
    for input in inputs:
        dirname = os.path.dirname(input)
        in_strs = input.split('/')[0].split('_')
        snake_content = 'conn = \'{}\'\nctsp = {}\nstp = {}\ng = {}\n'.format(in_strs[0], in_strs[1], in_strs[2], in_strs[3])
        snake_content += snake_template
        with open(os.path.join(dirname, 'Snakefile'), 'w') as f:
            f.write(snake_content)
