import time

def print_text(f, text):
    text = time.strftime("%Y-%m-%d %H:%M:%S ") + text.strip()
    print(text)
    f.write(text + "\n")