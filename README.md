# Grav

Количественный анализ элементарных клеточных автоматов (ECA): вычисление метрик, кластеризация правил и отдельные инструменты для поиска translating fronts.

## Quick Start

Через `make`:

```bash
make metrics
make cluster
```

Напрямую через Python:

```bash
python3 -m pip install -r requirements.txt
python3 -m grav.metrics --out-dir results
python3 -m grav.cluster --csv results/eca_metrics.csv --out-dir results
```

## Что делает проект

Для всех 256 правил ECA проект оценивает:

- среднюю энтропию состояния,
- чувствительность к одно-битовому возмущению,
- Lyapunov-like эвристику,
- грубый индикатор glider-like поведения.

После этого правила кластеризуются в пространстве признаков, а качество кластеризации проверяется через silhouette, устойчивость по `ARI` и пермутационный тест.

## Минимальные зависимости

```bash
python3 -m pip install -r requirements.txt
```

Для разработки и ноутбуков:

```bash
python3 -m pip install -r requirements-dev.txt
```

## Воспроизводимый запуск

Основной путь теперь идёт через пакет `grav`, а результаты по умолчанию сохраняются в `results/`.

```bash
python3 -m grav.metrics --runs 20 --N 200 --T 200 --out-dir results
python3 -m grav.cluster --csv results/eca_metrics.csv --out-dir results
```

То же самое через `Makefile`:

```bash
make reproduce
```

Полезные цели:

- `make metrics`
- `make cluster`
- `make test`

## Структура

- `grav/eca.py` — базовые ECA-примитивы и битовая реализация шага.
- `grav/nonlinear_core.py` — основная реализация вычисления метрик.
- `grav/cluster_core.py` — основная реализация кластеризации.
- `grav/metrics.py` — CLI для генерации метрик и базовых иллюстраций.
- `grav/cluster.py` — CLI для кластеризации и пермутационного теста.
- `grav/fronts.py` — пакетная точка входа для поиска translating fronts.
- `grav/research/` — одноразовые исследовательские скрипты, вынесенные из корня.
- `tests/` — минимальные тесты корректности.
- `docs/` — статический GitHub Pages сайт с обзором проекта, preprint-страницей и заметкой по Rule 69.

Совместимые обёртки в корне сохранены для старого CLI (`python3 nonlinear.py`, `python3 find_translating_fronts.py` и т.д.), но рекомендуемый запуск теперь идёт через `python -m grav...`.

## GitHub Pages

Сайт репозитория хранится в `docs/` и деплоится через GitHub Actions workflow `.github/workflows/pages.yml`.
Для публикации в настройках репозитория нужно выбрать GitHub Pages source = `GitHub Actions`.

## Артефакты

PNG, CSV, `fronts_out/` и `__pycache__/` считаются генерируемыми файлами и исключены из git через `.gitignore`.

## Проверка корректности

В репозиторий добавлен тест на эквивалентность двух реализаций шага ECA:

- посимвольной на кольце,
- битовой на целочисленном представлении.

Запуск:

```bash
python3 -m unittest discover -s tests -p 'test_*.py'
```

## Ограничения

- текущие метрики остаются эвристическими;
- glider/front detection не является исчерпывающей классификацией;
- большая часть исследовательских скриптов всё ещё живёт в корне и может быть позже перенесена в пакет без изменения логики.
